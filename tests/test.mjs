// node --test tests/
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import * as P from '../core/physics.mjs';
import * as PL from '../core/planner.mjs';

const D = JSON.parse(fs.readFileSync(new URL('../data/core_state.json', import.meta.url)));
const geom = P.buildGeometry();
const th = PL.defaultThimbles(geom);
const mkCore = () => ({ channels: D.channels.map(c => ({ ...c, bundles: c.bundles.slice(), xs: P.channelXS(c.bundles) })) });
const spec = { ...PL.DEFAULT_SPEC, rhoNominal: D.rhoNominal };

function homogeneousGeom(N, h) {
  const type = new Int8Array(N * N).fill(2), chanIndex = new Int32Array(N * N);
  const channels = [];
  for (let k = 0; k < N * N; k++) { chanIndex[k] = k; channels.push({ node: k, i: k % N, j: (k / N) | 0 }); }
  return { N, ctr: (N - 1) / 2, type, chanIndex, channels, h };
}

test('lattice: k-inf curve has a plutonium peak and crosses 1 near 1.5 n/kb; xs reproduce k-inf exactly', () => {
  assert.ok(P.kinf(0.15) > P.kinf(0) && P.kinf(0.15) > P.kinf(0.5));
  assert.ok(P.kinf(1.4) > 1 && P.kinf(1.6) < 1);
  const c = P.CONST;
  for (const w of [0, 0.3, 1, 1.8]) {
    const k = (c.nuSf1 + P.nuSf2(w) * c.S12 / P.Sa2(w)) / (c.Sa1 + c.S12);
    assert.ok(Math.abs(k - P.kinf(w)) < 1e-12);
  }
});

test('diffusion: bare homogeneous square converges to the analytic two-group eigenvalue under mesh refinement', () => {
  const c = P.CONST, w = 0.6, a = 400; // cm side
  const errs = [];
  for (const N of [10, 20, 40]) {
    const g = homogeneousGeom(N, a / N);
    const core = { channels: g.channels.map(() => ({ bundles: new Array(12).fill(w) })) };
    const m = P.materials(g, core, null, null, null);
    const k = P.solveDiffusion(g, m, null, { tol: 1e-10, maxOuter: 4000, inner: 6 }).k;
    // analytic with extrapolated side a + 2*2.13*D (group-averaged extrapolation)
    const ae1 = a + 2 * 2.13 * c.D1, ae2 = a + 2 * 2.13 * c.D2;
    const B1 = 2 * (Math.PI / ae1) ** 2 + m.Bz2, B2 = 2 * (Math.PI / ae2) ** 2 + m.Bz2;
    const kA = (c.nuSf1 + P.nuSf2(w) * c.S12 / (c.D2 * B2 + P.Sa2(w))) / (c.D1 * B1 + c.Sa1 + c.S12);
    errs.push(Math.abs(k - kA) / kA);
  }
  assert.ok(errs[2] < errs[0], `error should shrink with refinement: ${errs}`);
  assert.ok(errs[2] < 2e-3, `fine-mesh error ${errs[2]}`);
});

test('diffusion: symmetric uniform core gives a symmetric power map and powers sum to rated power', () => {
  const core = { channels: geom.channels.map(() => ({ bundles: new Array(12).fill(0.9) })) };
  const m = P.materials(geom, core, null, null, null);
  const sol = P.solveDiffusion(geom, m, null, { tol: 1e-10 });
  const pd = P.powerDistribution(geom, core, m, sol);
  const total = pd.chanMW.reduce((a, x) => a + x, 0);
  assert.ok(Math.abs(total - P.CONST.powerMWth) < 1e-6);
  const N = geom.N;
  geom.channels.forEach((ch, idx) => {
    const mirror = geom.chanIndex[ch.j * N + (N - 1 - ch.i)];
    assert.ok(Math.abs(pd.chanMW[idx] - pd.chanMW[mirror]) < 1e-5);
  });
  pd.bundleKW.forEach((b, i) => assert.ok(Math.abs(b.reduce((a, x) => a + x, 0) / 1000 - pd.chanMW[i]) < 1e-9));
});

test('refuelling: 8-bundle shift keeps the right bundles for both fuelling directions', () => {
  const core = { channels: [{ bundles: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11].map(x => x / 10), dir: 1 }, { bundles: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11].map(x => x / 10), dir: -1 }] };
  const d1 = P.refuel(core, 0, 5);
  assert.deepEqual(core.channels[0].bundles, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.2, 0.3]);
  assert.deepEqual(d1, [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1]);
  const d2 = P.refuel(core, 1, 5);
  assert.deepEqual(core.channels[1].bundles, [0.8, 0.9, 1, 1.1, 0, 0, 0, 0, 0, 0, 0, 0]);
  assert.deepEqual(d2, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]);
  assert.equal(core.channels[0].dir, 1); assert.equal(core.channels[0].lastRefuel, 5);
});

test('bateman: analytic chain matches RK4 integration and the burnout-limited saturation activity', () => {
  const iso = P.ISO.lu, N10 = 3.3e22, phi = 2e14, T = 9 * 86400;
  const r = P.bateman(iso, N10, 0, phi, T);
  const s1 = iso.targetSigma * 1e-24 * phi, lam = Math.LN2 / (iso.halfLifeDays * 86400), s2 = lam + iso.productSigma * 1e-24 * phi;
  let y = [N10, 0]; const dt = 600;
  const f = ([a, b]) => [-s1 * a, s1 * a - s2 * b];
  for (let t = 0; t < T; t += dt) {
    const k1 = f(y), k2 = f(y.map((v, i) => v + dt / 2 * k1[i])), k3 = f(y.map((v, i) => v + dt / 2 * k2[i])), k4 = f(y.map((v, i) => v + dt * k3[i]));
    y = y.map((v, i) => v + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
  }
  assert.ok(Math.abs(r.N2 - y[1]) / y[1] < 1e-8);
  const sat = P.bateman(iso, N10, 0, phi, 400 * 86400);
  const asym = s1 * N10 * Math.exp(-s1 * 400 * 86400) * lam / s2 / 3.7e10; // quasi-saturation with slow target burn
  assert.ok(Math.abs(sat.activityCi - asym) / asym < 1e-3);
  const hist = P.batemanHistory(iso, N10, new Array(9).fill(phi));
  assert.ok(Math.abs(hist.N2 - r.N2) / r.N2 < 1e-10);
});

test('Lu-177 cycle: chosen irradiation time is the brute-force optimum of net value per occupied day', () => {
  for (const phi of [1.2e14, 1.8e14, 2.4e14]) {
    const best = PL.bestLuCycle(phi, spec);
    for (let t = 1; t <= 28; t++) {
      const v = (PL.luCycle(phi, spec, t) * spec.prices.luPerCi - spec.prices.luCapsule) / (t + spec.handlingDays);
      assert.ok(v <= best.rate + 1e-9);
    }
  }
  assert.ok(PL.bestLuCycle(2.4e14, spec).t <= PL.bestLuCycle(1.0e14, spec).t, 'higher flux should not lengthen the optimal cycle');
});

test('absorbers: cobalt stringers cost reactivity and neighbouring stringers shadow each other (worth is sub-additive)', () => {
  const core = mkCore();
  const empty = th.map(() => 'empty');
  const r0 = PL.evaluate(geom, core, th, empty, spec).rho;
  const one = (i) => { const l = empty.slice(); l[i] = 'co'; return PL.evaluate(geom, core, th, l, spec).rho - r0; };
  const wa = one(0), wb = one(4); // T01 and T05 are neighbours (inner ring)
  const both = empty.slice(); both[0] = 'co'; both[4] = 'co';
  const wab = PL.evaluate(geom, core, th, both, spec).rho - r0;
  assert.ok(wa < 0 && wb < 0);
  assert.ok(Math.abs(wab) < Math.abs(wa + wb), `pair ${wab} vs sum ${wa + wb}`);
  const lu = empty.slice(); lu[0] = 'lu';
  assert.ok(Math.abs(PL.evaluate(geom, core, th, lu, spec).rho - r0) < 0.05 * Math.abs(wa), 'Yb capsule worth is small');
});

test('fuel sensitivity: parasitic absorption lowers discharge irradiation and raises annual fuel cost', () => {
  const s = PL.fuelSensitivity(geom, D.exitW0 / 2);
  assert.ok(s.drhoDw < 0 && s.exitPerMk > 0);
  const a = PL.annualFuel(spec, D.exitW0), b = PL.annualFuel(spec, D.exitW0 - 0.1);
  assert.ok(b.cost > a.cost);
});

test('optimiser: matches exhaustive search on a 6-thimble instance and respects capacity limits', () => {
  const core = mkCore();
  const th6 = th.slice(0, 6);
  const sp = { ...spec, luCapacity: 2, coInventory: 2 };
  const ctx = PL.makeContext(geom, core, th6, sp, D.refMW, D.exitW0);
  let bestEx = null;
  const kinds = ['empty', 'lu', 'co'];
  for (let code = 0; code < 3 ** 6; code++) {
    const l = []; let c = code; for (let i = 0; i < 6; i++) { l.push(kinds[c % 3]); c = Math.floor(c / 3); }
    if (l.filter(x => x === 'lu').length > 2 || l.filter(x => x === 'co').length > 2) continue;
    const s = PL.scoreLoading(geom, core, th6, l, sp, ctx);
    if (!bestEx || s.net > bestEx.net) bestEx = s;
  }
  const opt = PL.optimise(geom, core, th6, sp, ctx, { iters: 300, seed: 3 });
  // re-score both winners with a tight solver so warm-start tolerance noise cannot decide the comparison
  const tight = { ...sp, solverTol: 1e-10 };
  const tctx = PL.makeContext(geom, core, th6, tight, D.refMW, D.exitW0);
  const a = PL.scoreLoading(geom, core, th6, opt.best.loading, tight, tctx).net;
  const b = PL.scoreLoading(geom, core, th6, bestEx.loading, tight, tctx).net;
  assert.ok(a >= b * (1 - 1e-3), `${a} vs exhaustive ${b}`);
  assert.ok(opt.best.loading.filter(x => x === 'lu').length <= 2 && opt.best.loading.filter(x => x === 'co').length <= 2);
});

test('optimiser: full 16-thimble plan beats the highest-flux rule of thumb on the static objective', () => {
  const core = mkCore();
  const ctx = PL.makeContext(geom, core, th, spec, D.refMW, D.exitW0);
  const naive = PL.scoreLoading(geom, core, th, PL.naiveLoading(geom, core, th, spec, ctx), spec, ctx);
  const opt = PL.optimise(geom, core, th, spec, ctx, { iters: 500, seed: 7 });
  assert.ok(opt.best.net >= naive.net);
  assert.ok(opt.best.net >= opt.greedy.net);
  assert.ok(opt.best.feasible);
});

test('cobalt sweep: warm-started levels give non-decreasing net value as inventory grows', () => {
  const core = mkCore();
  const ctx = PL.makeContext(geom, core, th, spec, D.refMW, D.exitW0);
  let last = -Infinity, n = 0;
  for (const st of PL.sweepCobalt(geom, core, th, spec, ctx, { levels: [0, 3, 6, 10], iters: 150 })) {
    assert.ok(st.best.net >= last - 1e-3 * Math.abs(last), `level ${st.nCo}: ${st.best.net} < ${last}`); // within solver-tolerance noise
    assert.ok(st.best.loading.filter(x => x === 'co').length <= st.nCo);
    last = st.best.net; n++;
  }
  assert.equal(n, 4);
});

test('simulation: fuel-only operation holds zone-control fill in band and is deterministic', () => {
  const core = mkCore();
  const sp = { ...spec };
  const mk = () => ({ ...PL.makeContext(geom, core, th, sp, D.refMW, D.exitW0), exitTarget: D.exitTarget, dailyBurnMk: D.dailyBurnMk });
  const empty = th.map(() => 'empty');
  const a = PL.simulate(geom, core, th, empty, sp, mk(), { days: 20, warmup: 5 });
  const b = PL.simulate(geom, core, th, empty, sp, mk(), { days: 20, warmup: 5 });
  assert.equal(a.viol.lzc, 0);
  assert.ok(a.refuels > 10 && a.refuels < 45, `refuels ${a.refuels}`);
  assert.equal(a.refuels, b.refuels);
  assert.deepEqual(a.log.map(x => x.picks.join()), b.log.map(x => x.picks.join()));
});

test('data: committed equilibrium core reproduces CANDU-like discharge burnup and refuelling rate', () => {
  const bu = D.exitW0 * P.BURNUP_PER_NKB;
  assert.ok(bu > 6.5 && bu < 8.5, `burnup ${bu}`);
  assert.ok(D.refuelsPerDay > 0.8 && D.refuelsPerDay < 1.6);
  assert.equal(D.channels.length, geom.channels.length);
});
