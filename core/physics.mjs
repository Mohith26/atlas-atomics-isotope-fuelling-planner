// Core physics for the isotope-aware fuelling planner.
// Two-group, 2D (x-y) finite-difference neutron diffusion for a pressure-tube
// heavy-water reactor (one mesh node per lattice cell, axial leakage by buckling),
// bundle-level irradiation tracking, on-power refuelling, and isotope Bateman chains.
// Pure functions + plain objects so it runs unchanged in Node (tests) and the browser.

export const CONST = {
  pitch: 28.575,            // cm, lattice pitch (CANDU-style square lattice)
  length: 594.36,           // cm, active channel length (12 x 49.53 cm bundles)
  nBundles: 12,
  shift: 8,                 // bundles inserted per refuelling (8-bundle shift)
  axialExtrapolated: 650,   // cm, effective axial height incl. end reflector savings
  coreRadius: 7.95,         // in pitches -> number of fuel channels (see buildGeometry)
  reflectorThickness: 2.3,  // in pitches (~66 cm D2O)
  powerMWth: 1050,          // ~350 MWe at ~33% net efficiency (Utah MOU: ~350 MWe commercial unit)
  kappa: 3.204e-11,         // J per fission (200 MeV)
  nu: 2.5,                  // neutrons per fission (lumped)
  fuelFluxRatio: 0.45,      // fuel-to-cell thermal flux ratio; calibrated so bundle residence (~290 FPD to
                            // ~1.8 n/kb) and burn-down rate (~0.5 mk/FPD) match CANDU-6 experience
  kgUPerBundle: 19.2,
  // Fast group (fuel cells), held constant with irradiation
  D1: 1.2627, Sa1: 0.0010, S12: 0.0068, nuSf1: 0.0011,
  D2: 0.9384,
  // Heavy-water reflector
  refl: { D1: 1.31, D2: 0.85, Sa1: 0.0, S12: 0.0115, Sa2: 0.00012 },
};

// Lattice k-infinity vs bundle irradiation (n/kb) for natural-uranium 37-element-style
// fuel with equilibrium xenon: plutonium peak near 0.1-0.25 n/kb, k-inf = 1 near 1.5 n/kb.
// Shape follows published CANDU-6 lattice behaviour; values are representative, not a
// licensed lattice-code output.
export function kinf(w) { return 1.11 + 0.20 * w * Math.exp(-w / 0.25) - 0.075 * w; }
export function Sa2(w) { return 0.0040 * (1 + 0.04 * w); }
export function nuSf2(w) {
  const c = CONST;
  return (kinf(w) * (c.Sa1 + c.S12) - c.nuSf1) * Sa2(w) / c.S12;
}
// Irradiation (n/kb) to burnup (MWd/kgU): ~4,200 MWd/tU per n/kb for natural-U CANDU fuel
export const BURNUP_PER_NKB = 4.2; // MWd/kgU per n/kb

// Axial shape of the 12 bundles (chopped cosine over the extrapolated height), mean = 1
export function axialShape() {
  const c = CONST, n = c.nBundles, L = c.length, H = c.axialExtrapolated;
  const s = [];
  for (let b = 0; b < n; b++) {
    const z = (b + 0.5) * L / n - L / 2;
    s.push(Math.cos(Math.PI * z / H));
  }
  const m = s.reduce((a, x) => a + x, 0) / n;
  return s.map(x => x / m);
}
const AX = axialShape();
const AX2 = AX.map(x => x * x);
const AX2SUM = AX2.reduce((a, x) => a + x, 0);

// Square grid; channels where the node centre lies within coreRadius of the core centre,
// D2O reflector out to coreRadius + reflectorThickness, vacuum beyond.
export function buildGeometry(opts = {}) {
  const c = { ...CONST, ...opts };
  const R = c.coreRadius, Rr = R + c.reflectorThickness;
  const N = 2 * Math.ceil(Rr) + 2;
  const ctr = (N - 1) / 2;
  const type = new Int8Array(N * N); // 0 vacuum, 1 reflector, 2 fuel
  const chanIndex = new Int32Array(N * N).fill(-1);
  const channels = [];
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    const r = Math.hypot(i - ctr, j - ctr);
    const k = j * N + i;
    if (r <= R) {
      type[k] = 2; chanIndex[k] = channels.length;
      const row = String.fromCharCode(65 + j - Math.round(ctr - R));
      channels.push({ node: k, i, j, r, name: `${row}${String(i - Math.round(ctr - R) + 1).padStart(2, '0')}` });
    } else if (r <= Rr) type[k] = 1;
  }
  return { N, ctr, type, chanIndex, channels, h: c.pitch };
}

// Importance-weighted channel cross sections from its 12 bundle irradiations
export function channelXS(bundles) {
  let sa2 = 0, nsf2 = 0;
  for (let b = 0; b < bundles.length; b++) {
    sa2 += AX2[b] * Sa2(bundles[b]);
    nsf2 += AX2[b] * nuSf2(bundles[b]);
  }
  return { Sa2: sa2 / AX2SUM, nuSf2: nsf2 / AX2SUM };
}

// Thimble absorber: N*sigma (cm^2 per cm of height, smeared over the cell volume)
export function thimbleSigma(kind, spec) {
  if (kind === 'co') return spec.coAtoms * spec.coSigma * 1e-24 / (CONST.pitch * CONST.pitch * CONST.length);
  if (kind === 'lu') return spec.ybAtoms * spec.ybSigma * 1e-24 / (CONST.pitch * CONST.pitch * CONST.length);
  return 0;
}

// Assemble per-node material arrays for the current core + thimble loading
export function materials(geom, core, thimbles, loading, spec) {
  const c = CONST, n = geom.N * geom.N;
  const Bz2 = (Math.PI / c.axialExtrapolated) ** 2;
  const m = {
    D1: new Float64Array(n), D2: new Float64Array(n), R1: new Float64Array(n), S12: new Float64Array(n),
    A2: new Float64Array(n), F1: new Float64Array(n), F2: new Float64Array(n), Bz2,
  };
  for (let k = 0; k < n; k++) {
    const t = geom.type[k];
    if (t === 1) {
      const r = c.refl;
      m.D1[k] = r.D1; m.D2[k] = r.D2; m.S12[k] = r.S12; m.R1[k] = r.Sa1 + r.S12; m.A2[k] = r.Sa2;
    } else if (t === 2) {
      const ch = core.channels[geom.chanIndex[k]];
      const xs = ch.xs || channelXS(ch.bundles);
      m.D1[k] = c.D1; m.D2[k] = c.D2; m.S12[k] = c.S12; m.R1[k] = c.Sa1 + c.S12;
      m.A2[k] = xs.Sa2; m.F1[k] = c.nuSf1; m.F2[k] = xs.nuSf2;
    }
  }
  if (thimbles && loading) {
    thimbles.forEach((t, idx) => {
      const kind = loading[idx];
      if (kind && kind !== 'empty') m.A2[t.node] += thimbleSigma(kind, spec);
    });
  }
  return m;
}

// Two-group power iteration with Gauss-Seidel/SOR inner sweeps.
// Returns k-eff and unnormalised fluxes. Warm start through `guess`.
export function solveDiffusion(geom, m, guess = null, opts = {}) {
  const { N, type, h } = geom;
  const n = N * N, tol = opts.tol ?? 1e-7, maxOuter = opts.maxOuter ?? 400, inner = opts.inner ?? 3, sor = opts.sor ?? 1.5;
  const phi1 = guess ? Float64Array.from(guess.phi1) : new Float64Array(n);
  const phi2 = guess ? Float64Array.from(guess.phi2) : new Float64Array(n);
  if (!guess) for (let k = 0; k < n; k++) if (type[k]) { phi1[k] = 1; phi2[k] = 0.3; }
  // coupling coefficients (per unit volume)
  const nb = [-1, 1, -N, N];
  const c1 = new Float64Array(n * 4), c2 = new Float64Array(n * 4);
  const diag1 = new Float64Array(n), diag2 = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    if (!type[k]) continue;
    const i = k % N, j = (k / N) | 0;
    let d1 = m.R1[k] + m.D1[k] * m.Bz2, d2 = m.A2[k] + m.D2[k] * m.Bz2;
    for (let q = 0; q < 4; q++) {
      const ii = i + (q === 0 ? -1 : q === 1 ? 1 : 0), jj = j + (q === 2 ? -1 : q === 3 ? 1 : 0);
      const kk = k + nb[q];
      const inside = ii >= 0 && ii < N && jj >= 0 && jj < N && type[kk];
      if (inside) {
        const a1 = 2 * m.D1[k] * m.D1[kk] / (m.D1[k] + m.D1[kk]) / (h * h);
        const a2 = 2 * m.D2[k] * m.D2[kk] / (m.D2[k] + m.D2[kk]) / (h * h);
        c1[k * 4 + q] = a1; c2[k * 4 + q] = a2; d1 += a1; d2 += a2;
      } else {
        // zero flux at extrapolated distance 2.13 D beyond the outer face
        d1 += m.D1[k] / (h * (h / 2 + 2.13 * m.D1[k]));
        d2 += m.D2[k] / (h * (h / 2 + 2.13 * m.D2[k]));
      }
    }
    diag1[k] = d1; diag2[k] = d2;
  }
  let k_eff = guess?.k ?? 1.0;
  const src = new Float64Array(n), prev = new Float64Array(n);
  let fissOld = 0;
  for (let k = 0; k < n; k++) { src[k] = m.F1[k] * phi1[k] + m.F2[k] * phi2[k]; fissOld += src[k]; }
  for (let k = 0; k < n; k++) prev[k] = src[k] / fissOld;
  let it = 0, err = 1;
  for (it = 0; it < maxOuter; it++) {
    for (let s = 0; s < inner; s++) {
      for (let k = 0; k < n; k++) {
        if (!type[k]) continue;
        let r = src[k] / k_eff;
        for (let q = 0; q < 4; q++) { const a = c1[k * 4 + q]; if (a) r += a * phi1[k + nb[q]]; }
        phi1[k] += sor * (r / diag1[k] - phi1[k]);
        if (phi1[k] < 0) phi1[k] = 0;
      }
    }
    for (let s = 0; s < inner; s++) {
      for (let k = 0; k < n; k++) {
        if (!type[k]) continue;
        let r = m.S12[k] * phi1[k];
        for (let q = 0; q < 4; q++) { const a = c2[k * 4 + q]; if (a) r += a * phi2[k + nb[q]]; }
        phi2[k] += sor * (r / diag2[k] - phi2[k]);
        if (phi2[k] < 0) phi2[k] = 0;
      }
    }
    let fiss = 0, dmax = 0, smax = 0;
    for (let k = 0; k < n; k++) {
      const v = m.F1[k] * phi1[k] + m.F2[k] * phi2[k];
      fiss += v; src[k] = v;
    }
    const kNew = k_eff * fiss / fissOld;
    // source convergence: compare normalised shapes between outer iterations
    for (let k = 0; k < n; k++) { const a = src[k] / fiss, b = prev[k]; if (Math.abs(a - b) > dmax) dmax = Math.abs(a - b); if (a > smax) smax = a; prev[k] = a; }
    err = Math.max(Math.abs(kNew - k_eff), (opts.shapeFactor ?? 1) * dmax / (smax || 1));
    k_eff = kNew;
    // renormalise
    const sc = 1 / fiss * 1000;
    for (let k = 0; k < n; k++) { phi1[k] *= sc; phi2[k] *= sc; src[k] *= sc; }
    fissOld = 1000;
    if (err < tol && it > 5) break;
  }
  return { k: k_eff, phi1, phi2, iterations: it + 1, converged: err < tol };
}

// Normalise fluxes to reactor power; returns channel powers (MW), bundle powers (kW),
// and absolute fluxes (n/cm2/s)
export function powerDistribution(geom, core, m, sol) {
  const c = CONST, V = c.pitch * c.pitch * c.length;
  let total = 0;
  const raw = new Float64Array(geom.channels.length);
  geom.channels.forEach((ch, idx) => {
    const k = ch.node;
    const fr = (m.F1[k] * sol.phi1[k] + m.F2[k] * sol.phi2[k]) / c.nu; // fission rate density (relative)
    raw[idx] = fr * V; total += raw[idx];
  });
  const scale = c.powerMWth * 1e6 / (total * c.kappa); // multiply relative flux by this to get n/cm2/s
  const chanMW = Array.from(raw, x => x / total * c.powerMWth);
  const bundleKW = core.channels.map((ch, idx) => {
    let s = 0; const w = ch.bundles.map((om, b) => { const v = AX[b] * nuSf2(om); s += v; return v; });
    return w.map(v => v / s * chanMW[idx] * 1000);
  });
  const phi2abs = new Float64Array(sol.phi2.length), phi1abs = new Float64Array(sol.phi1.length);
  for (let k = 0; k < sol.phi2.length; k++) { phi2abs[k] = sol.phi2[k] * scale; phi1abs[k] = sol.phi1[k] * scale; }
  return { chanMW, bundleKW, phi1abs, phi2abs, scale };
}

// Advance bundle irradiations by `days` at the current flux
export function burn(geom, core, phi2abs, days) {
  const c = CONST, sec = days * 86400;
  geom.channels.forEach((g, idx) => {
    const ch = core.channels[idx];
    const f = phi2abs[g.node] * c.fuelFluxRatio * sec * 1e-21; // n/kb at mean axial position
    for (let b = 0; b < ch.bundles.length; b++) ch.bundles[b] += f * AX[b];
    ch.xs = channelXS(ch.bundles);
  });
}

// On-power refuel with an n-bundle shift. Each channel is always fuelled from the same end;
// adjacent channels are fuelled in opposite directions (checkerboard), as in CANDU cores.
export function refuel(core, idx, day) {
  const c = CONST, ch = core.channels[idx], n = c.nBundles, s = c.shift;
  const old = ch.bundles.slice();
  let discharged, nb;
  if (ch.dir >= 0) { nb = new Array(s).fill(0).concat(old.slice(0, n - s)); discharged = old.slice(n - s); }
  else { nb = old.slice(s).concat(new Array(s).fill(0)); discharged = old.slice(0, s); }
  ch.bundles = nb; ch.lastRefuel = day; ch.xs = channelXS(nb);
  return discharged;
}

export function cloneCore(core) {
  return { ...core, channels: core.channels.map(ch => ({ ...ch, bundles: ch.bundles.slice(), xs: ch.xs ? { ...ch.xs } : undefined })) };
}

// ---------------- Isotope chains ----------------
export const ISO = {
  lu: { // n.c.a. Lu-177 via 176Yb(n,g)177Yb (1.9 h) -> 177Lu; Yb-177 lumped into the decay
    targetSigma: 2.85, productSigma: 1000, halfLifeDays: 6.647, label: 'Lu-177 (n.c.a.)',
  },
  co: { targetSigma: 37.18, productSigma: 2.0, halfLifeDays: 5.2714 * 365.25, label: 'Co-60' },
};
const BQ_PER_CI = 3.7e10;

// Bateman for target -> product with product burnout under constant flux phi (n/cm2/s)
// over t seconds starting from N1(0)=N10, N2(0)=N20. Returns { N1, N2, activityCi }.
export function bateman(iso, N10, N20, phi, t) {
  const s1 = iso.targetSigma * 1e-24 * phi;
  const lam = Math.LN2 / (iso.halfLifeDays * 86400);
  const s2 = lam + iso.productSigma * 1e-24 * phi;
  const N1 = N10 * Math.exp(-s1 * t);
  let N2;
  if (Math.abs(s2 - s1) < 1e-30) N2 = N20 * Math.exp(-s2 * t) + s1 * N10 * t * Math.exp(-s1 * t);
  else N2 = N20 * Math.exp(-s2 * t) + s1 * N10 / (s2 - s1) * (Math.exp(-s1 * t) - Math.exp(-s2 * t));
  return { N1, N2, activityCi: lam * N2 / BQ_PER_CI };
}

// Integrate a piecewise-constant daily flux history through the chain
export function batemanHistory(iso, N10, fluxByDay) {
  let N1 = N10, N2 = 0;
  for (const phi of fluxByDay) { const r = bateman(iso, N1, N2, phi, 86400); N1 = r.N1; N2 = r.N2; }
  const lam = Math.LN2 / (iso.halfLifeDays * 86400);
  return { N1, N2, activityCi: lam * N2 / BQ_PER_CI };
}
