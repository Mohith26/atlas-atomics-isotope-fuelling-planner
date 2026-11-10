// Isotope-aware fuelling planner.
//  1. thimble assignment: simulated annealing over {empty, Lu-177 capsule, Co-59 stringer} per
//     isotope thimble, each candidate evaluated with a full two-group diffusion eigenvalue solve
//     (captures absorber worth, flux shadowing between neighbouring targets, and power peaking);
//  2. per-thimble Lu-177 cycle length that maximises activity per day of thimble occupancy;
//  3. daily channel selection by pre-simulation (each candidate refuel is solved before it is
//     chosen) against a reference power map rebased for the chosen target loading;
//  4. a horizon simulation that verifies the plan against operating limits and tallies realised
//     isotope yields and fuel consumption.
import * as P from './physics.mjs';

export const DEFAULT_SPEC = {
  coAtoms: 1500 / 58.933 * 6.022e23,      // 1.5 kg Co-59 per stringer
  coSigma: P.ISO.co.targetSigma,
  ybAtoms: 9.7 / 175.94 * 6.022e23,       // 10 g Yb, 97% Yb-176 per capsule
  ybSigma: P.ISO.lu.targetSigma,
  luCapacity: 6,                          // capsules per cycle (hot-cell throughput)
  coInventory: 10,                        // stringers on hand
  handlingDays: 2,                        // thimble dead time per Lu-177 cycle (unload, transfer, reload)
  prices: { luPerCi: 250, coPerCi: 3, fuelPerKgU: 250, luCapsule: 4000 }, // US$, illustrative assumptions
  capacityFactor: 0.9,
  limits: { maxChannelMW: 7.3, maxBundleKW: 935, cppf: 1.10, lzcLow: 0.20, lzcHigh: 0.80 }, // CANDU-6-style values; cppf is
                                          // replaced at run time by the fuel-only envelope (see app)
  lzcWorthMk: 7.0,                        // liquid zone control worth (mk) for fill-level conversion
  luCycleDays: 7,
  peakMarginMW: 0.05,                     // allowed static rise in peak channel power from targets
  minDwell: 60, candidates: 10,
};

// 16 vertical isotope thimbles at lattice interstitials (two rings); node coordinates on the grid
export function defaultThimbles(geom) {
  const pts = [[9, 9], [14, 9], [9, 14], [14, 14], [8, 11], [15, 12], [12, 8], [11, 15],
    [7, 7], [16, 7], [7, 16], [16, 16], [5, 12], [18, 11], [11, 5], [12, 18]];
  return pts.map(([i, j], idx) => {
    const node = j * geom.N + i;
    if (geom.type[node] !== 2) throw new Error('thimble outside fuelled region');
    return { id: `T${String(idx + 1).padStart(2, '0')}`, i, j, node, r: Math.hypot(i - geom.ctr, j - geom.ctr), channel: geom.chanIndex[node] };
  });
}

export const rho = k => (k - 1) / k * 1000; // mk

// Zones for power-shape control: inner/outer x four quadrants
export function zoneOf(geom, ch) {
  const q = (ch.i < geom.ctr ? 0 : 1) + (ch.j < geom.ctr ? 0 : 2);
  return (ch.r < 4.5 ? 0 : 4) + q;
}

export function evaluate(geom, core, thimbles, loading, spec, guess) {
  const m = P.materials(geom, core, thimbles, loading, spec);
  const sol = P.solveDiffusion(geom, m, guess, { tol: spec.solverTol ?? 1e-6, shapeFactor: 0.01 });
  const pd = P.powerDistribution(geom, core, m, sol);
  return { m, sol, pd, k: sol.k, rho: rho(sol.k) };
}

// Fuel-cost sensitivity: equilibrium exit-irradiation loss per mk of parasitic absorption,
// estimated numerically from uniform-irradiation cores around the operating average.
export function fuelSensitivity(geom, wAvg) {
  const kAt = w => {
    const core = { channels: geom.channels.map(() => ({ bundles: new Array(12).fill(w), dir: 1 })) };
    return P.solveDiffusion(geom, P.materials(geom, core, null, null, null), null, { tol: 1e-8 }).k;
  };
  const dw = 0.05;
  const drhoDw = (rho(kAt(wAvg + dw)) - rho(kAt(wAvg - dw))) / (2 * dw); // mk per n/kb of core-average
  return { drhoDw, exitPerMk: -2 / drhoDw }; // continuous fuelling: exit irradiation ~ 2x core average
}

export function annualFuel(spec, exitW) {
  const burnup = exitW * P.BURNUP_PER_NKB; // MWd/kgU
  const kgU = P.CONST.powerMWth * 365.25 * spec.capacityFactor / burnup;
  return { burnup, kgU, cost: kgU * spec.prices.fuelPerKgU };
}

// Lu-177 cycle: activity at end of irradiation for t days at flux phi; optimal t for throughput
export function luCycle(phi, spec, tDays) {
  const r = P.bateman(P.ISO.lu, spec.ybAtoms, 0, phi, tDays * 86400);
  return r.activityCi;
}
// Maximise net value per day of thimble occupancy: (Ci * price - capsule processing cost) / (t + handling)
export function bestLuCycle(phi, spec) {
  let best = null;
  for (let t = 1; t <= 28; t += 1) {
    const ci = luCycle(phi, spec, t), rate = (ci * spec.prices.luPerCi - spec.prices.luCapsule) / (t + spec.handlingDays);
    if (!best || rate > best.rate) best = { t, rate, ci };
  }
  return best;
}

// Static value of a thimble loading on the current core ($/yr), using diffusion fluxes
export function scoreLoading(geom, core, thimbles, loading, spec, ctx) {
  const ev = evaluate(geom, core, thimbles, loading, spec, ctx.guess);
  ctx.guess = ev.sol;
  const dRho = ev.rho - ctx.rho0; // negative for absorbers
  const days = 365.25 * spec.capacityFactor;
  let luCiYr = 0, coCiYr = 0, luCapsules = 0;
  const per = thimbles.map((t, idx) => {
    const phi = ev.pd.phi2abs[t.node];
    const kind = loading[idx];
    let ciYr = 0, cycle = null;
    if (kind === 'lu') {
      cycle = spec.optimiseCycle === false ? { t: spec.luCycleDays, ci: luCycle(phi, spec, spec.luCycleDays) } : bestLuCycle(phi, spec);
      const cycles = days / (cycle.t + spec.handlingDays);
      ciYr = cycle.ci * cycles; luCiYr += ciYr; luCapsules += cycles;
    } else if (kind === 'co') {
      // first-year accrual (Co-60 half-life 5.27 y), self-shielding folded into the cell flux
      const r = P.bateman(P.ISO.co, spec.coAtoms, 0, phi, days * 86400);
      ciYr = r.activityCi; coCiYr += ciYr;
    }
    return { id: t.id, kind, phi, ciYr, cycle };
  });
  const exitLoss = -dRho * ctx.sens.exitPerMk; // n/kb (positive = lower discharge irradiation)
  const exitW = ctx.exitW0 - exitLoss;
  const fuel = annualFuel(spec, exitW);
  const fuelPenalty = fuel.cost - ctx.fuel0.cost;
  const value = luCiYr * spec.prices.luPerCi + coCiYr * spec.prices.coPerCi - luCapsules * spec.prices.luCapsule;
  const maxMW = Math.max(...ev.pd.chanMW);
  const maxB = Math.max(...ev.pd.bundleKW.map(b => Math.max(...b)));
  // design check against the no-target core: a loading may not raise peak channel / bundle power
  // by more than the allowed margin (keeps the campaign inside the existing operating envelope)
  const mwCap = ctx.peakMW0 + spec.peakMarginMW, bCap = ctx.peakB0 + spec.peakMarginMW * 1000 / 12 * 1.6;
  const overLimit = Math.max(0, maxMW - mwCap) + Math.max(0, maxB - bCap) / 100;
  return { loading: loading.slice(), dRho, exitW, fuelPenalty, value, net: value - fuelPenalty - overLimit * 1e8,
    luCiYr, coCiYr, luCapsules, maxMW, maxB, feasible: overLimit === 0, per, k: ev.k, chanMW: ev.pd.chanMW };
}

// Deterministic PRNG (mulberry32)
export function rng(seed) {
  let a = seed >>> 0;
  return () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}

export function makeContext(geom, core, thimbles, spec, refMW, exitW0) {
  const ev0 = evaluate(geom, core, thimbles, thimbles.map(() => 'empty'), spec, null);
  const sens = fuelSensitivity(geom, exitW0 / 2);
  const peakMW0 = Math.max(...ev0.pd.chanMW), peakB0 = Math.max(...ev0.pd.bundleKW.map(b => Math.max(...b)));
  return { rho0: ev0.rho, guess: ev0.sol, sens, exitW0, fuel0: annualFuel(spec, exitW0), refMW, ev0, peakMW0, peakB0 };
}

const feasible = (loading, spec) =>
  loading.filter(x => x === 'lu').length <= spec.luCapacity && loading.filter(x => x === 'co').length <= spec.coInventory;

// Conventional rule: highest-flux thimbles get cobalt (inventory permitting), next ones Lu-177
export function naiveLoading(geom, core, thimbles, spec, ctx) {
  const order = thimbles.map((t, i) => [i, ctx.ev0.pd.phi2abs[t.node]]).sort((a, b) => b[1] - a[1]).map(x => x[0]);
  const loading = thimbles.map(() => 'empty');
  order.slice(0, spec.coInventory).forEach(i => loading[i] = 'co');
  order.slice(spec.coInventory, spec.coInventory + spec.luCapacity).forEach(i => loading[i] = 'lu');
  return loading;
}

// Greedy constructive start: repeatedly add the (thimble, isotope) with the best marginal net value
export function greedy(geom, core, thimbles, spec, ctx) {
  let cur = thimbles.map(() => 'empty');
  let curS = scoreLoading(geom, core, thimbles, cur, spec, ctx), evals = 1;
  for (;;) {
    let best = null;
    for (let i = 0; i < cur.length; i++) {
      if (cur[i] !== 'empty') continue;
      for (const kind of ['lu', 'co']) {
        const cand = cur.slice(); cand[i] = kind;
        if (!feasible(cand, spec)) continue;
        const sc = scoreLoading(geom, core, thimbles, cand, spec, ctx); evals++;
        if (!best || sc.net > best.net) best = sc;
      }
    }
    if (!best || best.net <= curS.net) break;
    cur = best.loading; curS = best;
  }
  return { best: curS, evals };
}

export function anneal(geom, core, thimbles, spec, ctx, opts = {}) {
  const R = rng(opts.seed ?? 7), iters = opts.iters ?? 900;
  const kinds = ['empty', 'lu', 'co'];
  let cur = opts.start ? opts.start.slice() : thimbles.map(() => 'empty');
  let curS = scoreLoading(geom, core, thimbles, cur, spec, ctx);
  let best = curS, evals = 1;
  const T0 = opts.T0 ?? Math.max(2e4, Math.abs(curS.net) * 0.01);
  const trace = [];
  for (let it = 0; it < iters; it++) {
    const T = T0 * Math.pow(0.002, it / iters);
    const cand = cur.slice();
    if (R() < 0.5) { // swap two thimbles' contents
      const a = Math.floor(R() * cand.length), b = Math.floor(R() * cand.length);
      [cand[a], cand[b]] = [cand[b], cand[a]];
    } else {
      const a = Math.floor(R() * cand.length);
      cand[a] = kinds[(kinds.indexOf(cand[a]) + 1 + Math.floor(R() * 2)) % 3];
    }
    if (!feasible(cand, spec) || cand.join() === cur.join()) continue;
    const s = scoreLoading(geom, core, thimbles, cand, spec, ctx); evals++;
    if (s.net > curS.net || R() < Math.exp((s.net - curS.net) / T)) { cur = cand; curS = s; }
    if (s.net > best.net) best = s;
    if (it % 30 === 0) trace.push({ it, net: curS.net, best: best.net });
  }
  return { best, evals, trace };
}

// Best-improvement local search over single changes and pairwise swaps until no move helps
export function polish(geom, core, thimbles, spec, ctx, start) {
  let cur = start, evals = 0;
  const kinds = ['empty', 'lu', 'co'];
  for (let round = 0; round < 50; round++) {
    let best = cur;
    const n = cur.loading.length;
    const tryL = (l) => { if (!feasible(l, spec)) return; const s = scoreLoading(geom, core, thimbles, l, spec, ctx); evals++; if (s.net > best.net + 1e-6) best = s; };
    for (let i = 0; i < n; i++) for (const k of kinds) if (k !== cur.loading[i]) { const l = cur.loading.slice(); l[i] = k; tryL(l); }
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (cur.loading[i] !== cur.loading[j]) { const l = cur.loading.slice(); [l[i], l[j]] = [l[j], l[i]]; tryL(l); }
    if (best === cur) break;
    cur = best;
  }
  return { best: cur, evals };
}

// Full optimiser: greedy construction (or a supplied warm start), annealing, then local-search
// polish. Any extra feasible starts (e.g. the optimum for a smaller cobalt inventory) are polished
// too, so the result is never worse than a known feasible plan.
export function optimise(geom, core, thimbles, spec, ctx, opts = {}) {
  let evals = 0;
  const g = opts.warm ? null : greedy(geom, core, thimbles, spec, ctx);
  if (g) evals += g.evals;
  const start = opts.warm ? scoreLoading(geom, core, thimbles, opts.warm, spec, ctx) : g.best;
  const a = anneal(geom, core, thimbles, spec, ctx, { ...opts, start: start.loading });
  evals += a.evals;
  let best = a.best.net >= start.net ? a.best : start;
  const pol = polish(geom, core, thimbles, spec, ctx, best); evals += pol.evals; best = pol.best;
  return { best, greedy: g ? g.best : start, evals, trace: a.trace };
}

// Sweep the cobalt inventory 0..max, warm-starting each level from the previous optimum
// (feasible at the larger inventory), so net value is non-decreasing along the sweep.
export function* sweepCobalt(geom, core, thimbles, spec, ctx, opts = {}) {
  const levels = opts.levels ?? [0, 2, 4, 6, 8, 10].filter(n => n < spec.coInventory).concat([spec.coInventory]);
  let prev = null, total = 0;
  const pts = [];
  for (const nCo of levels) {
    const sp = { ...spec, coInventory: nCo };
    const r = optimise(geom, core, thimbles, sp, ctx, { iters: opts.iters ?? 320, seed: 11 + nCo, warm: prev ? prev.loading : null });
    let best = r.best;
    if (prev) { const re = scoreLoading(geom, core, thimbles, prev.loading, sp, ctx); if (re.net > best.net) best = re; }
    total += r.evals; prev = best;
    pts.push({ nCo, best });
    yield { nCo, best, pts, evals: total };
  }
}

// Pareto sweep: best loading for each Co-60 stringer count (0..inventory)
export function paretoByCobalt(geom, core, thimbles, spec, ctx, opts = {}) {
  const pts = [];
  for (let nCo = 0; nCo <= spec.coInventory; nCo += opts.step ?? 2) {
    const sp = { ...spec, coInventory: nCo };
    const start = thimbles.map(() => 'empty');
    const r = optimise(geom, core, thimbles, sp, ctx, { iters: opts.iters ?? 150, seed: 11 + nCo });
    pts.push({ nCo, ...r.best });
  }
  return pts;
}

// Reference power map for a target loading: time-average reference scaled by the static
// redistribution the targets cause (so selection does not fight a permanent, intended shape change)
export function rebasedReference(geom, core, thimbles, loading, spec, ctx) {
  const e = evaluate(geom, core, thimbles, loading, spec, null);
  return ctx.refMW.map((r, i) => r * e.pd.chanMW[i] / ctx.ev0.pd.chanMW[i]);
}

// ---------------- daily channel selection + horizon simulation ----------------
export function lzcFill(rhoMk, spec) { return 0.5 + (rhoMk - spec.rhoNominal) / spec.lzcWorthMk; }

export function chooseRefuels(geom, core, thimbles, loading, spec, day, ctx) {
  const picks = [];
  let ev = evaluate(geom, core, thimbles, loading, spec, ctx.guess); ctx.guess = ev.sol;
  const target = spec.rhoNominal; // aim for 50% LZC fill at end of day
  const burnDrop = ctx.dailyBurnMk ?? 0.35;
  let guard = 0;
  while (ev.rho - burnDrop < target - 0.6 && guard++ < 6) {
    const zMW = new Float64Array(8), zRef = new Float64Array(8);
    // candidate list: eligible channels ranked by fuelling-engineer rules
    const cands = core.channels.map((ch, i) => {
      const g = geom.channels[i];
      const wbar = ch.bundles.reduce((a, x) => a + x, 0) / ch.bundles.length;
      const since = day - ch.lastRefuel;
      const powRatio = ev.pd.chanMW[i] / ctx.refMW[i];
      return { i, prio: wbar / ctx.exitTarget[i] - 0.5 * (powRatio - 1), since, wbar };
    }).filter(c => c.since >= spec.minDwell && !picks.includes(c.i))
      .sort((a, b) => b.prio - a.prio).slice(0, spec.candidates ?? 10);
    let best = null;
    for (const c of cands) {
      const trial = P.cloneCore(core); P.refuel(trial, c.i, day);
      const e = evaluate(geom, trial, thimbles, loading, spec, ev.sol);
      const ratios = e.pd.chanMW.map((p, i) => p / ctx.refMW[i]);
      const cppf = Math.max(...ratios), maxMW = Math.max(...e.pd.chanMW);
      const maxB = Math.max(...e.pd.bundleKW.map(b => Math.max(...b)));
      zMW.fill(0); zRef.fill(0);
      geom.channels.forEach((g, i) => { const z = zoneOf(geom, g); zMW[z] += e.pd.chanMW[i]; zRef[z] += ctx.refMW[i]; });
      let tilt = 0; for (let z = 0; z < 8; z++) tilt += (zMW[z] / zRef[z] - 1) ** 2;
      let J = 40 * tilt + 30 * Math.max(0, cppf - 1.04) ** 2 * 100 - 0.4 * c.prio;
      J += 1e3 * (Math.max(0, cppf - spec.limits.cppf) + Math.max(0, maxMW - spec.limits.maxChannelMW) + Math.max(0, maxB - spec.limits.maxBundleKW) / 100);
      if (!best || J < best.J) best = { ...c, J, e, cppf, maxMW, maxB };
    }
    if (!best) break;
    P.refuel(core, best.i, day); picks.push(best.i);
    ev = best.e; ctx.guess = ev.sol;
  }
  return { picks, ev };
}

// Simulate `days` of operation (after `warmup` unscored days). Generator form yields progress so a
// browser can stay responsive; simulate() runs it to completion.
export function simulate(geom, core0, thimbles, loading, spec, ctx, opts) {
  const it = simulateIter(geom, core0, thimbles, loading, spec, ctx, opts);
  for (;;) { const r = it.next(); if (r.done) return r.value; }
}
export function* simulateIter(geom, core0, thimbles, loading, spec, ctx, opts) {
  const core = P.cloneCore(core0);
  const days = opts.days, warm = opts.warmup ?? 0;
  const log = [], discharged = [];
  const luIdx = loading.map((k, i) => k === 'lu' ? i : -1).filter(i => i >= 0);
  const cycleLen = {};
  luIdx.forEach(i => cycleLen[i] = Math.max(1, Math.round(opts.cycleDays?.[i] ?? spec.luCycleDays)));
  const cap = {}; // per-capsule state
  luIdx.forEach(i => cap[i] = { start: 0, flux: [] });
  let luCi = 0, luHarvests = 0, coN2 = loading.map(() => 0), coN1 = loading.map(k => k === 'co' ? spec.coAtoms : 0);
  let refuels = 0, viol = { lzc: 0, cppf: 0, mw: 0, bundle: 0 };
  let maxCPPF = 0, maxMW = 0, maxB = 0, lzcMin = 9, lzcMax = -9;
  for (let d = -warm; d < days; d++) {
    const counting = d >= 0;
    const { picks, ev } = chooseRefuels(geom, core, thimbles, loading, spec, d, ctx);
    // burn one day at the post-refuel flux
    P.burn(geom, core, ev.pd.phi2abs, 1);
    if (counting) {
      refuels += picks.length;
      const ratios = ev.pd.chanMW.map((p, i) => p / ctx.refMW[i]);
      const cppf = Math.max(...ratios), mw = Math.max(...ev.pd.chanMW);
      const mb = Math.max(...ev.pd.bundleKW.map(b => Math.max(...b)));
      const fill = lzcFill(ev.rho, spec);
      maxCPPF = Math.max(maxCPPF, cppf); maxMW = Math.max(maxMW, mw); maxB = Math.max(maxB, mb);
      lzcMin = Math.min(lzcMin, fill); lzcMax = Math.max(lzcMax, fill);
      if (fill < spec.limits.lzcLow || fill > spec.limits.lzcHigh) viol.lzc++;
      if (cppf > spec.limits.cppf) viol.cppf++;
      if (mw > spec.limits.maxChannelMW) viol.mw++;
      if (mb > spec.limits.maxBundleKW) viol.bundle++;
      // isotopes: Lu-177 capsules irradiate for cycleLen days, then the thimble is idle for handlingDays
      for (const i of luIdx) {
        const c = cap[i];
        if (c.skip > 0) { c.skip--; continue; }
        c.flux.push(ev.pd.phi2abs[thimbles[i].node]);
        if (c.flux.length >= cycleLen[i]) {
          const r = P.batemanHistory(P.ISO.lu, spec.ybAtoms, c.flux);
          luCi += r.activityCi; luHarvests++;
          c.flux = []; c.skip = Math.round(spec.handlingDays);
        }
      }
      loading.forEach((k, i) => {
        if (k !== 'co') return;
        const r = P.bateman(P.ISO.co, coN1[i], coN2[i], ev.pd.phi2abs[thimbles[i].node], 86400);
        coN1[i] = r.N1; coN2[i] = r.N2;
      });
      yield { day: d, of: days };
      log.push({ day: d, picks: picks.map(i => geom.channels[i].name), pickIdx: picks.slice(), rho: ev.rho, fill, cppf, maxMW: mw, maxBundleKW: mb, chanMW: d === days - 1 ? ev.pd.chanMW : null });
    }
  }
  const coCi = coN2.reduce((a, n) => a + Math.LN2 / (P.ISO.co.halfLifeDays * 86400) * n / 3.7e10, 0);
  return { log, refuels, bundles: refuels * P.CONST.shift, luCi, luHarvests, coCi, viol, maxCPPF, maxMW, maxB, lzcMin, lzcMax, finalCore: core };
}
