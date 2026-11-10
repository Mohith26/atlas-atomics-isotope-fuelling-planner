import * as P from './core/physics.mjs';
import * as PL from './core/planner.mjs';

// match the site's viewport-driven scale variables (same formula as atlasatomics.com)
(function syncScale() {
  const r = document.documentElement, w = innerWidth, h = innerHeight;
  const k = Math.min(1.6, Math.max(0.34, Math.min(w / 1920, h / 1080))), vmax = Math.max(w, h);
  r.style.setProperty('--k', k.toFixed(4));
  r.style.setProperty('--halo-scale', ((1.6 * 1.7 * vmax) / 800).toFixed(4));
  r.style.setProperty('--blob-scale', ((1.2 * 1.7 * vmax) / 800).toFixed(4));
})();

const $ = (s) => document.querySelector(s);
const tick = () => new Promise(r => setTimeout(r, 0));
const fmt = (x, d = 0) => x.toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
const money = (x, d = 1) => (x < 0 ? '-' : '') + '$' + fmt(Math.abs(x) / 1e6, d) + 'M';
const signed = (x, d = 1, unit = '') => (x >= 0 ? '+' : '-') + fmt(Math.abs(x), d) + unit;
const DAYS = 42, WARM = 30;

const S = { view: 'plan' };

async function boot() {
  const D = await (await fetch('data/core_state.json')).json();
  S.D = D;
  S.geom = P.buildGeometry();
  S.th = PL.defaultThimbles(S.geom);
  document.querySelectorAll('[data-k="bu"]').forEach(e => e.textContent = fmt(D.exitW0 * P.BURNUP_PER_NKB, 2));
  document.querySelectorAll('[data-k="rf"]').forEach(e => e.textContent = fmt(D.refuelsPerDay, 2));
  document.querySelectorAll('[data-k="warm"]').forEach(e => e.textContent = WARM);
  for (const [i, o] of [['#in-co', '#out-co'], ['#in-lu', '#out-lu']]) $(i).addEventListener('input', () => $(o).textContent = $(i).value);
  $('#controls').addEventListener('submit', (e) => { e.preventDefault(); run(); });
  document.querySelectorAll('#view-seg button').forEach(b => b.addEventListener('click', () => {
    document.querySelectorAll('#view-seg button').forEach(x => x.classList.toggle('on', x === b));
    S.view = b.dataset.v; drawMap();
  }));
  await run();
}

function freshCore() {
  return { channels: S.D.channels.map(c => ({ ...c, bundles: c.bundles.slice(), xs: P.channelXS(c.bundles) })) };
}

function readSpec() {
  return {
    ...PL.DEFAULT_SPEC, rhoNominal: S.D.rhoNominal,
    coInventory: +$('#in-co').value, luCapacity: +$('#in-lu').value,
    prices: { ...PL.DEFAULT_SPEC.prices, luPerCi: +$('#in-plu').value, coPerCi: +$('#in-pco').value, fuelPerKgU: +$('#in-pf').value },
  };
}

function status(t, done = false) { $('#chip-status').textContent = t; $('.chip.live').classList.toggle('done', done); }

async function run() {
  const token = (S.token = (S.token || 0) + 1);
  $('#run').disabled = true; status('Solving');
  await tick();
  const D = S.D, geom = S.geom, th = S.th;
  const spec = readSpec();
  const core = freshCore();
  const ctx = { ...PL.makeContext(geom, core, th, spec, D.refMW, D.exitW0), exitTarget: D.exitTarget, dailyBurnMk: D.dailyBurnMk };
  S.spec = spec; S.core = core; S.ctx = ctx;
  $('#chip-k').textContent = fmt(ctx.ev0.k, 5);
  $('#chip-rho').textContent = fmt(ctx.rho0, 2) + ' mk';
  $('#chip-bu').textContent = fmt(D.exitW0 * P.BURNUP_PER_NKB, 2) + ' MWd/kgU';
  $('#chip-rf').textContent = fmt(D.refuelsPerDay, 2);

  const ruleL = PL.naiveLoading(geom, core, th, spec, ctx);
  S.rule = PL.scoreLoading(geom, core, th, ruleL, { ...spec, optimiseCycle: false }, ctx);
  await tick();
  // plan = top of a warm-started sweep over cobalt inventory (also gives the trade-off curve)
  const pts = [];
  for (const st of PL.sweepCobalt(geom, core, th, spec, ctx)) {
    pts.push({ nCo: st.nCo, s: st.best });
    S.plan = st.best; S.pts = pts;
    $('#thimble-note').textContent = `${fmt(st.evals)} diffusion solves`;
    status(`Planning ${pts.length}`);
    renderKPIs(); drawMap(); renderThimbles(); drawPareto(pts);
    await tick(); if (token !== S.token) return;
  }
  paretoText(pts);
  status('Verifying');
  await verify(token);
  if (token !== S.token) return;
  status('Plan ready', true);
  $('#run').disabled = false;
}

function renderKPIs() {
  const p = S.plan, r = S.rule, base = S.ctx;
  const netP = p.value - p.fuelPenalty, netR = r.value - r.fuelPenalty;
  const cards = [
    ['Net isotope value / yr', money(netP, 2), `<b>${signed((netP - netR) / 1e6, 2, 'M')}</b> vs rule of thumb (${money(netR, 2)})`],
    ['Lu-177 at EOI / yr', fmt(p.luCiYr / 1000, 1) + '<small>kCi</small>', `<b>${signed((p.luCiYr / r.luCiYr - 1) * 100, 0, '%')}</b> · ${fmt(p.luCapsules, 0)} capsules / yr`],
    ['Co-60 accrued / yr', fmt(p.coCiYr / 1e6, 2) + '<small>MCi</small>', `<b class="neg">${signed((p.coCiYr / Math.max(r.coCiYr, 1) - 1) * 100, 0, '%')}</b> · ${p.loading.filter(x => x === 'co').length} of ${S.spec.coInventory} stringers used`],
    ['Reactivity cost', fmt(-p.dRho, 2) + '<small>mk</small>', `fuel ${money(p.fuelPenalty, 2)}/yr · rule ${fmt(-r.dRho, 2)} mk, ${money(r.fuelPenalty, 2)}`],
    ['Peak channel power', fmt(p.maxMW, 2) + '<small>MW</small>', `target-free core ${fmt(base.peakMW0, 2)} MW · rule ${fmt(r.maxMW, 2)}`],
  ];
  $('#kpis').innerHTML = cards.map(([c, v, d]) => `<div class="kpi"><div class="caption">${c}</div><div class="v">${v}</div><div class="d">${d}</div></div>`).join('');
}

// ---------- core map ----------
const CELL = 25.5, C0 = 280;
const pos = (i, j) => [C0 + (i - S.geom.ctr) * CELL, C0 + (j - S.geom.ctr) * CELL];
function ramp(t) { // dark taupe -> pale -> copper
  t = Math.max(0, Math.min(1, t));
  const stops = [[0, [70, 58, 60, .55]], [0.55, [226, 214, 206, .9]], [1, [255, 118, 62, 1]]];
  let a = stops[0], b = stops[1];
  if (t > stops[1][0]) { a = stops[1]; b = stops[2]; }
  const u = (t - a[0]) / (b[0] - a[0]);
  const c = a[1].map((x, k) => x + (b[1][k] - x) * u);
  return `rgba(${c[0] | 0},${c[1] | 0},${c[2] | 0},${c[3].toFixed(2)})`;
}
function drawMap() {
  if (!S.plan) return;
  const sc = S.view === 'plan' ? S.plan : S.rule, geom = S.geom;
  const mw = sc.chanMW, lo = 2.5, hi = 7.4;
  const Rr = (P.CONST.coreRadius + P.CONST.reflectorThickness) * CELL;
  let s = `<circle cx="${C0}" cy="${C0}" r="${Rr}" fill="rgba(14,10,16,.18)" stroke="rgba(255,255,255,.15)" stroke-width="1"/>`;
  s += `<circle cx="${C0}" cy="${C0}" r="${(P.CONST.coreRadius + .5) * CELL}" fill="none" stroke="rgba(255,255,255,.15)" stroke-width="1" stroke-dasharray="2 4"/>`;
  for (const rr of [3.54, 6.4]) s += `<circle cx="${C0}" cy="${C0}" r="${rr * CELL}" fill="none" stroke="rgba(255,255,255,.08)" stroke-width="1"/>`;
  geom.channels.forEach((ch, idx) => {
    const [x, y] = pos(ch.i, ch.j);
    s += `<circle class="ch" data-c="${idx}" cx="${x}" cy="${y}" r="${CELL * 0.36}" fill="${ramp((mw[idx] - lo) / (hi - lo))}"/>`;
  });
  S.th.forEach((t, idx) => {
    const [x, y] = pos(t.i + 0.5, t.j + 0.5), kind = sc.loading[idx], d = 6.2;
    const fill = kind === 'lu' ? '#ff763e' : kind === 'co' ? '#dfe3ea' : 'rgba(14,10,16,.55)';
    const stroke = kind === 'empty' ? 'rgba(255,255,255,.6)' : 'rgba(14,10,16,.7)';
    s += `<path class="th" data-t="${idx}" d="M${x} ${y - d}L${x + d} ${y}L${x} ${y + d}L${x - d} ${y}Z" fill="${fill}" stroke="${stroke}" stroke-width="1"/>`;
  });
  // axis labels (rows A.. / columns 01..)
  const rows = [...new Set(geom.channels.map(c => c.j))].sort((a, b) => a - b), cols = [...new Set(geom.channels.map(c => c.i))].sort((a, b) => a - b);
  rows.forEach((j, k) => { const [, y] = pos(0, j); s += `<text x="16" y="${y + 3}" fill="rgba(255,255,255,.4)" font-size="7.5" letter-spacing="1">${String.fromCharCode(65 + k)}</text>`; });
  cols.forEach((i, k) => { if (k % 3) return; const [x] = pos(i, 0); s += `<text x="${x - 5}" y="14" fill="rgba(255,255,255,.4)" font-size="7.5">${String(k + 1).padStart(2, '0')}</text>`; });
  const svg = $('#map'); svg.innerHTML = s;
  $('#legend').innerHTML = `<span>Channel MW ${lo}<i class="ramp" style="background:linear-gradient(90deg,${ramp(0)},${ramp(.55)},${ramp(1)})"></i>${hi}</span>
    <span><i class="sw" style="background:#ff763e"></i>Lu-177 capsule</span><span><i class="sw" style="background:#dfe3ea"></i>Co-59 stringer</span><span><i class="sw" style="border:1px solid rgba(255,255,255,.6)"></i>Empty thimble</span>`;
  const tip = $('#tip');
  svg.onmousemove = (e) => {
    const t = e.target, wrap = svg.getBoundingClientRect();
    let html = null;
    if (t.classList.contains('ch')) {
      const i = +t.dataset.c, ch = S.core.channels[i], g = geom.channels[i];
      const wbar = ch.bundles.reduce((a, x) => a + x, 0) / 12;
      html = `<span class="caption">Channel ${g.name}</span>${fmt(mw[i], 2)} MW · ref ${fmt(S.ctx.refMW[i], 2)} MW<br>mean irradiation ${fmt(wbar, 3)} n/kb · last fuelled ${-ch.lastRefuel} FPD ago`;
    } else if (t.classList.contains('th')) {
      const i = +t.dataset.t, p = sc.per[i];
      html = `<span class="caption">Thimble ${p.id} · ${p.kind === 'lu' ? 'Lu-177' : p.kind === 'co' ? 'Co-60' : 'empty'}</span>thermal flux ${fmt(p.phi / 1e14, 3)}×10¹⁴ n/cm²s` + (p.kind !== 'empty' ? `<br>${fmt(p.ciYr, 0)} Ci / yr` + (p.cycle ? ` · ${p.cycle.t}-day cycle` : '') : '');
    }
    if (!html) { tip.style.opacity = 0; return; }
    tip.innerHTML = html; tip.style.opacity = 1;
    tip.style.left = Math.min(e.clientX - wrap.left + 14, wrap.width - 250) + 'px'; tip.style.top = (e.clientY - wrap.top + 14) + 'px';
  };
  svg.onmouseleave = () => tip.style.opacity = 0;
}

function tag(kind, changed) {
  const lab = kind === 'lu' ? 'Lu-177' : kind === 'co' ? 'Co-60' : 'Empty';
  return `<span class="tag ${kind}${changed ? ' changed' : ''}"><i></i>${lab}</span>`;
}
function renderThimbles() {
  const p = S.plan, r = S.rule;
  const rows = S.th.map((t, i) => ({ t, i, phi: S.ctx.ev0.pd.phi2abs[t.node] })).sort((a, b) => b.phi - a.phi);
  $('#thimbles').innerHTML = `<thead><tr><th>Thimble</th><th>Ring</th><th class="num">Flux ×10¹⁴</th><th>Rule of thumb</th><th>Planner</th><th class="num">Cycle</th><th class="num">Ci / yr</th></tr></thead><tbody>` +
    rows.map(({ t, i, phi }) => {
      const pk = p.loading[i], rk = r.loading[i], pp = p.per[i];
      return `<tr class="${pk !== rk ? 'hl' : ''}"><td class="id">${t.id}</td><td>${t.r < 5 ? 'Inner' : 'Outer'}</td><td class="num">${fmt(phi / 1e14, 2)}</td><td>${tag(rk)}</td><td>${tag(pk, pk !== rk)}</td><td class="num">${pp.cycle ? pp.cycle.t + ' d' : '-'}</td><td class="num">${pk === 'empty' ? '-' : fmt(pp.ciYr, 0)}</td></tr>`;
    }).join('') + '</tbody>';
}

// ---------- verification ----------
async function simulateAsync(loading, cycleDays, token, label) {
  const ctx = { ...S.ctx, refMW: PL.rebasedReference(S.geom, S.core, S.th, loading, S.spec, S.ctx), guess: null };
  const it = PL.simulateIter(S.geom, S.core, S.th, loading, S.spec, ctx, { days: DAYS, warmup: WARM, cycleDays });
  let n = 0;
  for (;;) {
    const r = it.next();
    if (r.done) return r.value;
    if (++n % 4 === 0) { $('#sim-progress').textContent = `${label} · FPD ${r.value.day + 1} / ${DAYS}`; await tick(); if (token !== S.token) return null; }
  }
}
async function verify(token) {
  const cyc = {}; S.plan.per.forEach((p, i) => { if (p.cycle) cyc[i] = p.cycle.t; });
  const empty = S.th.map(() => 'empty');
  const runs = [['Fuel only', empty, {}], ['Rule of thumb', S.rule.loading, {}], ['Planner', S.plan.loading, cyc]];
  S.sims = [];
  for (const [label, l, c] of runs) {
    const r = await simulateAsync(l, c, token, label);
    if (!r) return;
    S.sims.push({ label, ...r });
    renderVerify();
  }
  $('#sim-progress').textContent = `${WARM} + ${DAYS} FPD · ${S.sims.reduce((a, s) => a + s.refuels, 0)} refuels simulated`;
}
function realised(s) {
  const sp = S.spec, base = S.sims[0];
  const fuelKg = (s.refuels - base.refuels) * P.CONST.shift * P.CONST.kgUPerBundle;
  const v = s.luCi * sp.prices.luPerCi + s.coCi * sp.prices.coPerCi - s.luHarvests * sp.prices.luCapsule - fuelKg * sp.prices.fuelPerKgU;
  return v * 365.25 * sp.capacityFactor / DAYS;
}
function renderVerify() {
  const sims = S.sims;
  const best = sims.length === 3 ? (realised(sims[2]) >= realised(sims[1]) ? 2 : 1) : -1;
  $('#verify').innerHTML = `<thead><tr><th>Strategy</th><th class="num">Lu-177 Ci</th><th class="num">Co-60 kCi</th><th class="num">Refuels</th><th class="num">Peak MW</th><th class="num">LZC fill</th><th class="num">Net / yr</th></tr></thead><tbody>` +
    sims.map((s, i) => `<tr class="${i === best ? 'best' : ''}"><td>${s.label}</td><td class="num">${fmt(s.luCi, 0)}</td><td class="num">${fmt(s.coCi / 1000, 0)}</td><td class="num">${s.refuels}</td><td class="num">${fmt(s.maxMW, 2)}</td><td class="num">${fmt(s.lzcMin * 100, 0)}-${fmt(s.lzcMax * 100, 0)}%</td><td class="num">${i === 0 ? '-' : money(realised(s), 2)}</td></tr>`).join('') + '</tbody>';
  drawPeakChart();
}
function drawPeakChart() {
  const W = 620, H = 150, L = 40, R = 10, T = 12, B = 24, sims = S.sims;
  const all = sims.flatMap(s => s.log.map(d => d.maxMW));
  const lo = Math.floor((Math.min(...all) - 0.05) * 10) / 10, hi = Math.ceil((Math.max(...all) + 0.05) * 10) / 10;
  const x = d => L + (d / (DAYS - 1)) * (W - L - R), y = v => T + (1 - (v - lo) / (hi - lo)) * (H - T - B);
  let s = '';
  for (let k = 0; k <= 4; k++) { const v = lo + (hi - lo) * k / 4; s += `<line x1="${L}" x2="${W - R}" y1="${y(v)}" y2="${y(v)}" stroke="rgba(255,255,255,.1)"/><text x="${L - 6}" y="${y(v) + 3}" text-anchor="end" font-size="8" fill="rgba(255,255,255,.5)">${fmt(v, 2)}</text>`; }
  for (let d = 0; d < DAYS; d += 7) s += `<text x="${x(d)}" y="${H - 6}" text-anchor="middle" font-size="8" fill="rgba(255,255,255,.5)">${d}</text>`;
  const styles = [['rgba(255,255,255,.5)', '3 3'], ['#dfe3ea', ''], ['#ff763e', '']];
  sims.forEach((sm, i) => { s += `<path d="${sm.log.map((d, k) => (k ? 'L' : 'M') + x(k).toFixed(1) + ' ' + y(d.maxMW).toFixed(1)).join('')}" fill="none" stroke="${styles[i][0]}" stroke-width="${i === 2 ? 1.8 : 1.2}" stroke-dasharray="${styles[i][1]}"/>`; });
  s += `<text x="${L}" y="9" font-size="8" fill="rgba(255,255,255,.55)" letter-spacing="1.5">PEAK CHANNEL POWER (MW) BY FPD</text>`;
  $('#peak-chart').innerHTML = `<svg viewBox="0 0 ${W} ${H}">${s}</svg>`;
  $('#peak-legend').innerHTML = sims.map((sm, i) => `<span><i class="line" style="border-top:${i === 2 ? 2 : 1}px ${i === 0 ? 'dashed' : 'solid'} ${styles[i][0]}"></i>${sm.label}</span>`).join('');
}

// ---------- trade-off ----------
function paretoText(pts) {
  const maxCo = S.spec.coInventory, zero = pts[0].s, p = S.plan, r = S.rule;
  const ring = (l) => { let inner = 0, outer = 0; l.forEach((k, i) => { if (k === 'co') (S.th[i].r < 5 ? inner++ : outer++); }); return [inner, outer]; };
  const [pi, po] = ring(p.loading), [ri, ro] = ring(r.loading);
  const used = pi + po, netP = p.value - p.fuelPenalty, netR = r.value - r.fuelPenalty;
  $('#pareto-note').textContent = `${pts.length} plans · Co-59 inventory swept`;
  $('#pareto-text').innerHTML = `Lu-177 alone earns <b>${money(zero.value - zero.fuelPenalty, 2)}/yr</b>. With up to ${maxCo} stringers the planner loads <b>${used}</b> (${pi} inner, ${po} outer) for <b>${fmt(-p.dRho, 1)} mk</b>; the rule of thumb loads ${ri + ro} (${ri} inner, ${ro} outer) for ${fmt(-r.dRho, 1)} mk, pays <b>${money(r.fuelPenalty - p.fuelPenalty, 2)}/yr</b> more in fuel and nets <b>${money(netP - netR, 2)}/yr</b> less.`;
}
function drawPareto(pts) {
  const W = 560, H = 250, L = 46, R = 16, T = 14, B = 30;
  const all = [...pts.map(p => p.s), S.rule];
  const lastN = pts[pts.length - 1].nCo;
  const xs = all.map(s => s.fuelPenalty / 1e6), ys = all.map(s => s.value / 1e6);
  const x0 = 0, x1 = Math.max(...xs) * 1.12 + 0.05, y0 = Math.floor(Math.min(...ys) * 0.92), y1 = Math.ceil(Math.max(...ys) * 1.04);
  const x = v => L + (v - x0) / (x1 - x0) * (W - L - R), y = v => T + (1 - (v - y0) / (y1 - y0)) * (H - T - B);
  let s = '';
  for (let k = 0; k <= 4; k++) { const v = y0 + (y1 - y0) * k / 4; s += `<line x1="${L}" x2="${W - R}" y1="${y(v)}" y2="${y(v)}" stroke="rgba(255,255,255,.1)"/><text x="${L - 6}" y="${y(v) + 3}" text-anchor="end" font-size="8" fill="rgba(255,255,255,.5)">$${fmt(v, 0)}M</text>`; }
  for (let k = 0; k <= 4; k++) { const v = x0 + (x1 - x0) * k / 4; s += `<text x="${x(v)}" y="${H - 14}" text-anchor="middle" font-size="8" fill="rgba(255,255,255,.5)">$${fmt(v, 1)}M</text>`; }
  s += `<text x="${(L + W - R) / 2}" y="${H - 2}" text-anchor="middle" font-size="7.5" letter-spacing="1.5" fill="rgba(255,255,255,.5)">FUEL PENALTY / YR</text>`;
  s += `<text x="${L}" y="9" font-size="7.5" letter-spacing="1.5" fill="rgba(255,255,255,.5)">GROSS ISOTOPE VALUE / YR</text>`;
  // iso-net lines
  const sorted = pts.slice().sort((a, b) => a.s.fuelPenalty - b.s.fuelPenalty);
  s += `<path d="${sorted.map((p, k) => (k ? 'L' : 'M') + x(p.s.fuelPenalty / 1e6).toFixed(1) + ' ' + y(p.s.value / 1e6).toFixed(1)).join('')}" fill="none" stroke="#ff763e" stroke-width="1.4"/>`;
  // group levels that converged to the same plan so labels do not overprint
  const groups = [];
  pts.forEach(p => { const g = groups.find(g => g.key === p.s.loading.join()); if (g) g.n.push(p.nCo); else groups.push({ key: p.s.loading.join(), n: [p.nCo], s: p.s }); });
  groups.forEach(g => {
    const cx = x(g.s.fuelPenalty / 1e6), cy = y(g.s.value / 1e6), chosen = g.n.includes(lastN);
    const used = g.s.loading.filter(k => k === 'co').length;
    s += `<circle cx="${cx}" cy="${cy}" r="${chosen ? 5.5 : 3.6}" fill="${chosen ? '#ff763e' : 'rgba(14,10,16,.8)'}" stroke="#ff763e" stroke-width="1.3"/>`;
    s += `<text x="${cx + 9}" y="${cy + 13}" font-size="7.5" fill="rgba(255,255,255,.7)" letter-spacing=".5">${used} CO${g.n.length > 1 ? ` (CAP ${g.n[0]}-${g.n[g.n.length - 1]})` : ''}</text>`;
  });
  const rx = x(S.rule.fuelPenalty / 1e6), ry = y(S.rule.value / 1e6);
  s += `<path d="M${rx} ${ry - 6}L${rx + 6} ${ry}L${rx} ${ry + 6}L${rx - 6} ${ry}Z" fill="#dfe3ea"/><text x="${rx - 9}" y="${ry - 9}" text-anchor="end" font-size="7.5" fill="rgba(255,255,255,.8)">RULE OF THUMB</text>`;
  $('#pareto').innerHTML = `<svg viewBox="0 0 ${W} ${H}">${s}</svg>`;
}

boot().catch(err => { console.error(err); status('Error'); });
