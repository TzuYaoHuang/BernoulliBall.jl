# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Julia simulation studying the "Bernoulli ball" effect: a ball floating in a vertical jet, and
whether/how the jet's flow field exerts a lateral restoring force that keeps the ball centred. There
is no package source (no `src/`, no exported module) — the entire project is the simulation script(s)
under `examples/`, driven by `Project.toml`/`Manifest.toml`. Treat this as a research script repo, not
a library.

Two variants of the experiment live in `examples/`:
- `TwoD_BallOnJet.jl` — 2D, custom GLMakie/Pathlines viz (pressure + pathline "smoke" panels) plus a
  lateral force coefficient (`Cy`) plot.
- `ThreeD_BallOnJet.jl` — 3D, same jet/ball/`BiotSimulation` physics generalised to a square
  cross-section domain (circular nozzle, ball offset along one transverse axis only). Visualisation
  and time-stepping are both handled by WaterLily's `viz!` (from its GLMakie extension), rendering
  λ2-criterion vortex cores rather than pressure/pathlines — the same recipe as
  `ThreeD_TaylorGreenVortex.jl` in [WaterLily-Examples](https://github.com/WaterLily-jl/WaterLily-Examples).
  Because `viz!` owns the run loop, there's no manual `mom_step!` loop; the lateral force coefficient
  is instead sampled from inside `viz!`'s own per-frame callback (see walkthrough below).

## Commands

Run from the repo root (Julia 1.12.7, pinned in `Manifest.toml`):

```sh
# one-time / after editing Project.toml — resolves the [sources] git deps too
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# run the 2D simulation (renders ball_on_jet.mp4 and restoring_force.png into the repo root)
julia --project=. examples/TwoD_BallOnJet.jl

# run the 3D simulation (renders ball_on_jet_3d.mp4 into the repo root; needs a display for GLMakie)
julia --project=. examples/ThreeD_BallOnJet.jl
```

There is no test suite, linter, or CI config in this repo.

GPU: the script auto-detects CUDA (`CUDA.functional() ? CuArray : Array`) and falls back to CPU
`Array` automatically — nothing to configure to run on CPU-only machines.

## Architecture

The simulation is built entirely on the WaterLily.jl CFD ecosystem; understanding either script
requires knowing how these pieces compose (all pulled in via `[sources]` in `Project.toml`, i.e.
installed straight from their GitHub repos rather than the general registry). Note `Jet` and `ball()`
are each defined independently in the 2D and 3D scripts (not shared via a module) — a change to one
does not propagate to the other.

- **WaterLily** — the base immersed-boundary Navier-Stokes solver (`Simulation`, `mom_step!`, `AutoBody`).
- **BiotSavartBCs** — replaces WaterLily's default far-field boundary handling. `BiotSimulation` sets
  the domain's open boundaries (all faces except those listed in `nonbiotfaces`) from the interior
  vorticity field via a Biot-Savart integral, rather than a slip wall or simple convective exit. In
  both scripts' `ball()`, only the floor (`-1`, the jet inlet) is excluded from this treatment.
- **LilyPad** — supplies the 2nd-order departure-point particle advection scheme used by `Pathlines.update!`
  (2D script only).
- **Pathlines** (2D script only) — particle-tracer visualisation: `Particles` seeds tracers across the
  domain, `update!` advects them each step, and `PathlineCanvas`/`fade!`/`draw!` rasterise them into a
  fading, speed-coloured "dye" image — independent of and at a different resolution than the solver grid.
- **GLMakie** vs **Plots** (2D script) — both are loaded and both export overlapping verbs
  (`contourf!`, etc.), so Makie calls are always qualified (`GLMakie.contourf!`) to disambiguate;
  `Plots` is used only for the final force-vs-time line plot.
- **`viz!`** (3D script only) — from WaterLily's GLMakie package extension (`ext/WaterLilyMakieExt.jl`,
  loaded automatically once `GLMakie` is `using`'d). Owns both time-stepping (`sim_step!` internally)
  and rendering: each frame it calls the supplied `f(cpu_array, sim)` to fill a CPU buffer, then volume-
  renders it (`algorithm=:absorption` here) alongside the immersed body. This replaces the 2D script's
  manual `mom_step!`/`GLMakie.record` loop entirely — don't combine the two patterns in one script.
- **Adapt/CUDA** — `mem = CuArray` or `Array` is threaded through `ball()` to pick the simulation's
  array backend. The custom `Jet` boundary-condition closure needs a manual
  `Adapt.adapt_structure(to, j::Jet) = j` because its three fields share a single type parameter, which
  breaks `Adapt`'s generic per-field closure adaptor — this is a load-bearing workaround, not
  boilerplate, if you touch `Jet` in either script.

### `examples/TwoD_BallOnJet.jl` walkthrough

1. `Jet` — a custom `uBC(i,x,t)` inflow: uniform vertical velocity `U` inside radius `r` of the domain
   centreline, zero elsewhere and zero tangentially (a blower nozzle in an otherwise closed floor). Also
   initialises the jet column at `t=0`.
2. `ball(...)` — builds the `BiotSimulation`: a disk body (`AutoBody`, SDF-defined) of diameter `D`,
   offset `off·R` from the jet centreline, inside a domain of size `(H*D, W*D)`.
3. Sets up `Particles`/`PathlineCanvas` for the smoke visualisation, and precomputes the static body SDF
   (`WaterLily.measure_sdf!`) to mask pressure inside the (physically meaningless) immersed body.
4. Two-panel `GLMakie` figure: pressure field + body outline (left), pathline canvas + body outline
   (right). Both panels are `permutedims`'d because WaterLily's grid dim 1 is the flow (vertical)
   direction but Makie plots a matrix's dim 1 along x.
5. Main loop (inside `GLMakie.record` → `ball_on_jet.mp4`): advances the flow (`mom_step!`), advects
   particles, updates the pathline canvas, and accumulates the lateral force coefficient
   `Cy = -2·pressure_force(sim)[2]/L` over time.
6. Final `Plots.plot` of `Cy` vs. `tU/L` saved to `restoring_force.png`. Per the inline comment, `Cy`
   is sharply negative on startup (Bernoulli suction pulling the ball back toward the axis) then
   oscillates with the shed wake — its long-time mean is small and sign-sensitive to `off`/`Re`, not a
   settled restoring constant.

### `examples/ThreeD_BallOnJet.jl` walkthrough

Same physical setup as the 2D script, generalised to 3D: `Jet`'s inflow region becomes a circular
nozzle in the (dim2,dim3) plane (centred, since that cross-section is square), and the ball is a
sphere offset `off·R` along dim2 only, leaving dim3 centred on the jet axis. `nonbiotfaces=(-1,)` still
excludes just the floor. Default `D` is smaller than the 2D script's (2⁴ vs 2⁵) since a 3D grid at the
same resolution is an order of magnitude more cells.

Visualisation is a single `viz!(sim; f=λ₂!, duration, step, algorithm=:absorption, colormap=:Reds,
video="ball_on_jet_3d.mp4")` call — no manual loop. `λ₂!` computes `log10(-λ2)` into `sim.flow.σ` (the
λ2 vortex-core criterion of Jeong & Hussain 1995, scaled by `L/U`) and copies it to the CPU buffer
`viz!` renders; this is the same recipe used for `ThreeD_TaylorGreenVortex.jl` in WaterLily-Examples.
Unlike the 2D script there's no `sym` kwarg to `viz!` (mirror symmetry, used e.g. for a centred jet),
since the off-axis ball breaks domain symmetry.

`viz!` calls `f` once per rendered frame, after stepping to that frame's time — at that point the flow
field is fully projected (divergence-free), unlike mid-step during `udf` (which WaterLily calls on
unprojected predictor/corrector state, before `project!`). So `λ₂!` doubles as the force-sampling hook:
it also pushes `Cy = -2·pressure_force(sim)[2]/sim.L` and `sim_time(sim)` into module-level `Cy`/`t_F`
arrays each frame, and after `viz!` returns, `Plots.plot(t_F,Cy,...)` saves `restoring_force_3d.png` —
same idea as the 2D script's force plot, just sampled once per render `step` instead of every internal
solver step (`viz!` exposes no finer-grained hook).

Key physical parameters (all keyword args of `ball()` in both scripts): `D` grid resolution/ball
diameter, `Re` Reynolds number, `U` jet velocity, `off` ball's off-axis offset (fraction of radius),
`H`/`W` domain size in ball diameters, `h` ball's initial height, `r` jet radius (fraction of ball
radius).

Generated/local artifacts (`ball_on_jet.mp4`, `ball_on_jet_3d.mp4`, `restoring_force.png`,
`restoring_force_3d.png`, `Manifest.toml`) are gitignored — they exist locally but are never tracked,
so don't expect `git status` to show them.
