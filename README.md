# Isotope-aware fuelling planner (concept demo for Atlas Atomics)

**Live demo:** https://mohith26.github.io/atlas-atomics-isotope-fuelling-planner/

Unofficial concept demo by Mohith Gajjela. Not affiliated with Atlas Atomics. All core data is synthetic.

## What it is
A heavy-water reactor that makes power *and* isotopes has to decide, every cycle, what goes into each isotope position: a Lu-177 capsule (Yb-176 target, short half-life, high value, almost no reactivity cost), a Co-59 stringer (long-lived Co-60, a real absorber that costs reactivity and therefore fuel burnup), or nothing. Those decisions interact: cobalt shadows neighbouring capsules, raises local power elsewhere, and changes what the on-power refuelling has to do every day.

This page plans that jointly for a synthetic 208-channel, ~1,050 MWth pressure-tube core (sized to the ~350 MWe commercial unit described in Atlas's public Utah MOU), styled as a "Fuelling / Isotope planner" page in an Atlas-branded console.

## How it works
- **Core physics (`core/physics.mjs`)**: two-group 2D finite-difference diffusion (one node per 28.575 cm lattice cell, D2O reflector, axial buckling), power iteration with SOR inner sweeps; 12 bundles per channel with importance-weighted cross sections; 8-bundle-shift on-power refuelling; analytic Bateman chains for 176Yb(n,γ)177Lu with Lu-177 burnout and 59Co(n,γ)60Co.
- **Planner (`core/planner.mjs`)**:
  - thimble assignment: greedy construction, simulated annealing, then swap/flip local search; every candidate loading is scored with a full eigenvalue solve, so absorber worth, shadowing and peaking are priced in. Reactivity cost is converted to discharge-burnup loss through a numerically derived dρ/dω of the model core, then to fuel cost;
  - per-capsule Lu-177 irradiation time that maximises net value per day of thimble occupancy (handling dead time and capsule processing cost included);
  - a warm-started sweep over cobalt inventory (0..N stringers) that traces the isotope-value vs fuel-penalty trade-off;
  - verification: each strategy runs 30 unscored + 42 scored full-power days of daily channel selection by pre-simulation (each candidate refuel is solved before it is picked) against a reference power map rebased for the target loading; realised Lu-177/Co-60 yields, refuels, peak channel power and liquid-zone fill are tallied.
- **Data (`scripts/gen_core.mjs`)**: seeded generator that burns a randomised core in for 640 days with the same fuelling rules to reach an equilibrium-like state (exit burnup ~7.4 MWd/kgU, ~1.1 refuels per FPD, CANDU-6-like).

Everything on the page is computed in the browser from `data/core_state.json` by this code; nothing is hard-coded.

## Results on the default synthetic case (prices are illustrative inputs)
The rule of thumb (cobalt in the hottest thimbles, Lu-177 in the next, fixed 7-day cycles) is compared with the planner. On the default inputs the planner moves cobalt mostly to the outer ring, puts Lu-177 in the hottest inner thimbles with 6-day cycles, uses ~5 mk of reactivity instead of ~12.7 mk, and nets ~$1.4M/yr more on the static model; the 42-FPD simulation confirms the direction (more Lu-177, fewer refuels, peak channel power inside the fuel-only envelope). Exact figures are on the page and change with the inputs.

## Tests
`node --test tests/test.mjs` (13 tests): k-inf/xs consistency; bare-square eigenvalue converges to the analytic two-group value under mesh refinement; symmetric power map and power normalisation; refuelling shift for both directions; Bateman vs RK4 and saturation limit; Lu-177 cycle optimum vs brute force; cobalt worth is sub-additive (shadowing); fuel sensitivity sign; optimiser equals exhaustive search on a 6-thimble instance; full plan beats the rule of thumb; warm-started cobalt sweep is non-decreasing; simulation determinism and zone-control band; equilibrium-data sanity.

## Honest limitations
- Not Atlas data and not Atlas's design: Atlas is in stealth and has published nothing about its lattice, channels, isotope positions or fuelling scheme. I assumed a CANDU-lineage pressure-tube, on-power-refuelled core because that is the mainstream natural-uranium heavy-water design.
- Lattice constants follow a representative natural-uranium k-inf(irradiation) curve, not a licensed lattice-code (WIMS/DRAGON) output; one diffusion node per lattice cell; 2D with axial buckling; no xenon transients, adjusters or spatial zone control (zone fill is a uniform-worth proxy).
- Isotope targets are smeared into the cell's thermal absorption; no target self-shielding model. Cobalt value uses first-year Co-60 accrual.
- Prices (Lu-177 $/Ci, Co-60 $/Ci, fuel $/kgU, capsule processing) are illustrative inputs, not market data.
- The model core's refuelling ripple is larger than a real CANDU's, so power limits are compared against the fuel-only envelope rather than licensing limits.

## Novelty note
CANDU refuelling-channel selection optimisation is well studied (MILP, genetic algorithms, REFUEL/OPTIMA-style rules) and CANDU isotope production exists (Co-60 adjusters, Lu-177 at Bruce). I found no published tool or Atlas material that co-optimises isotope-target assignment, capsule timing and the fuelling plan against one core model; Atlas itself has shipped no software publicly.
