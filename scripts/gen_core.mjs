// Seeded generator for the synthetic equilibrium core used by the demo.
// Starts from random channel phases, then runs the fuel-only planner for a long burn-in so the
// core reaches an equilibrium-like irradiation distribution. The time-averaged channel powers
// over the last 150 days become the reference power distribution (as in CANDU fuel management).
// Usage: node scripts/gen_core.mjs  -> data/core_state.json
import fs from 'node:fs';
import * as P from '../core/physics.mjs';
import * as PL from '../core/planner.mjs';

const SEED = 20261109;
const R = PL.rng(SEED);
const geom = P.buildGeometry();
const thimbles = PL.defaultThimbles(geom);
const spec = { ...PL.DEFAULT_SPEC, minDwell: 60, candidates: 10, rhoNominal: 3.5 };
const Rc = P.CONST.coreRadius;

// design reference shape: flattened inner region, rolling off toward the edge
const shape = geom.channels.map(c => c.r <= 4.5 ? 1 : 1 - 0.45 * ((c.r - 4.5) / (Rc - 4.5)) ** 1.5);
const sum = shape.reduce((a, x) => a + x, 0);
let refMW = shape.map(s => s / sum * P.CONST.powerMWth);
const exitTarget = geom.channels.map(c => 1 + 0.15 * (1 - c.r / Rc));

const DAYS = 520, AVG = 150;
// random phases: bundles 0-7 inserted u*dwell ago, 8-11 one dwell earlier
const core = { channels: geom.channels.map((c, i) => {
  const u = R(), wd = 0.95 * (0.8 + 0.4 * shape[i]);
  const b = new Array(12).fill(0).map((_, k) => (k < 8 ? u : 1 + u) * wd * (0.85 + 0.3 * R()));
  return { bundles: b, dir: (geom.channels[i].i + geom.channels[i].j) % 2 ? 1 : -1, lastRefuel: -DAYS - Math.round(u * 190) };
}) };
core.channels.forEach(ch => ch.xs = P.channelXS(ch.bundles));

const loading = thimbles.map(() => 'empty');
const ctx = { refMW, exitTarget, guess: null, dailyBurnMk: 0.5 };
const acc = new Float64Array(geom.channels.length);
const disch = [];
let refuels = 0;
const t0 = Date.now();
for (let d = -DAYS; d < 0; d++) {
  const { picks, ev } = PL.chooseRefuels(geom, core, thimbles, loading, spec, d, ctx);
  if (d >= -AVG) { ev.pd.chanMW.forEach((p, i) => acc[i] += p / AVG); refuels += picks.length; }
  P.burn(geom, core, ev.pd.phi2abs, 1);
  if ((d + DAYS) % 40 === 0) console.log('day', d, 'rho', ev.rho.toFixed(2), 'refuels', picks.length, 'elig', core.channels.filter(c => d - c.lastRefuel >= spec.minDwell).length, ((Date.now() - t0) / 1000).toFixed(1) + 's');
}
refMW = Array.from(acc);
// measure equilibrium discharge irradiation & burn-down rate with the new reference
ctx.refMW = refMW;
let dis = [], rf = 0, rhoDrop = [];
for (let d = 0; d < 120; d++) {
  const before = P.cloneCore(core);
  const { picks, ev } = PL.chooseRefuels(geom, core, thimbles, loading, spec, d, ctx);
  picks.forEach(i => { const ch = before.channels[i]; const out = ch.dir >= 0 ? ch.bundles.slice(4) : ch.bundles.slice(0, 8); dis.push(...out); });
  rf += picks.length;
  const e0 = PL.evaluate(geom, core, thimbles, loading, spec, ctx.guess);
  P.burn(geom, core, e0.pd.phi2abs, 1);
  const e1 = PL.evaluate(geom, core, thimbles, loading, spec, e0.sol);
  rhoDrop.push(e0.rho - e1.rho); ctx.guess = e1.sol;
}
const exitW0 = dis.reduce((a, x) => a + x, 0) / dis.length;
const dailyBurnMk = rhoDrop.reduce((a, x) => a + x, 0) / rhoDrop.length;
const ev = PL.evaluate(geom, core, thimbles, loading, spec, null);
const out = {
  meta: { seed: SEED, generated: 'scripts/gen_core.mjs', burnInDays: DAYS + 120, note: 'Synthetic equilibrium core; not Atlas Atomics data.' },
  rhoNominal: spec.rhoNominal, exitW0, dailyBurnMk, refuelsPerDay: rf / 120,
  refMW, exitTarget, rhoNow: ev.rho,
  channels: core.channels.map(ch => ({ bundles: ch.bundles.map(x => +x.toFixed(5)), dir: ch.dir, lastRefuel: ch.lastRefuel - 120 })),
};
fs.writeFileSync(new URL('../data/core_state.json', import.meta.url), JSON.stringify(out));
console.log({ exitW0, burnup: exitW0 * P.BURNUP_PER_NKB, dailyBurnMk, refuelsPerDay: rf / 120, rhoNow: ev.rho, maxRef: Math.max(...refMW), minRef: Math.min(...refMW), secs: (Date.now() - t0) / 1000 });
