# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Julia simulation studying the "Bernoulli ball" effect: a ball floating in a vertical jet, and
whether/how the jet's flow field exerts a lateral restoring force that keeps the ball centred. There
is no package source (no `src/`, no exported module) — the entire project is the simulation script(s)
under `examples/`, driven by `Project.toml`/`Manifest.toml`. Treat this as a research script repo, not
a library.

## Commands

Run from the repo root (Julia 1.12.7, pinned in `Manifest.toml`):

```sh
# one-time / after editing Project.toml — resolves the [sources] git deps too
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# run the simulation (renders ball_on_jet.mp4 and restoring_force.png into the repo root)
julia --project=. examples/TwoD_BallOnJet.jl
```

There is no test suite, linter, or CI config in this repo.

GPU: the script auto-detects CUDA (`CUDA.functional() ? CuArray : Array`) and falls back to CPU
`Array` automatically — nothing to configure to run on CPU-only machines.

## Architecture

The simulation is built entirely on the WaterLily.jl CFD ecosystem; understanding the script requires
knowing how these pieces compose (all pulled in via `[sources]` in `Project.toml`, i.e. installed
straight from their GitHub repos rather than the general registry):

- **WaterLily** — the base immersed-boundary Navier-Stokes solver (`Simulation`, `mom_step!`, `AutoBody`).
- **BiotSavartBCs** — replaces WaterLily's default far-field boundary handling. `BiotSimulation` sets
  the domain's open boundaries (all faces except those listed in `nonbiotfaces`) from the interior
  vorticity field via a Biot-Savart integral, rather than a slip wall or simple convective exit. In
  `ball()`, only the floor (`-1`, the jet inlet) is excluded from this treatment.
- **LilyPad** — supplies the 2nd-order departure-point particle advection scheme used by `Pathlines.update!`.
- **Pathlines** — particle-tracer visualisation: `Particles` seeds tracers across the domain, `update!`
  advects them each step, and `PathlineCanvas`/`fade!`/`draw!` rasterise them into a fading,
  speed-coloured "dye" image — independent of and at a different resolution than the solver grid.
- **GLMakie** vs **Plots** — both are loaded and both export overlapping verbs (`contourf!`, etc.), so
  Makie calls are always qualified (`GLMakie.contourf!`) to disambiguate; `Plots` is used only for the
  final force-vs-time line plot.
- **Adapt/CUDA** — `mem = CuArray` or `Array` is threaded through `ball()` to pick the simulation's
  array backend. The custom `Jet` boundary-condition closure needs a manual
  `Adapt.adapt_structure(to, j::Jet) = j` because its three fields share a single type parameter, which
  breaks `Adapt`'s generic per-field closure adaptor — this is a load-bearing workaround, not
  boilerplate, if you touch `Jet`.

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

Key physical parameters (all keyword args of `ball()`): `D` grid resolution/ball diameter, `Re`
Reynolds number, `U` jet velocity, `off` ball's off-axis offset (fraction of radius), `H`/`W` domain
size in ball diameters, `h` ball's initial height, `r` jet radius (fraction of ball radius).

Generated/local artifacts (`ball_on_jet.mp4`, `restoring_force.png`, `Manifest.toml`) are gitignored —
they exist locally but are never tracked, so don't expect `git status` to show them.
