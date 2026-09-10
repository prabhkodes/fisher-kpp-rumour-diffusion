# fisher-kpp-rumour-diffusion

A rumour spreading across a 2-D population, solved as a Fisher-KPP reaction-diffusion problem in
parallel with PETSc. Awareness starts at one point and forms a travelling wave that moves outward at a
constant speed you can predict analytically and then check against the simulation.

**Stack:** C · PETSc (`DMDA`, `TS`) · MPI · ParaView / VisIt

**Equation:**

```
∂u/∂t = D ∇²u + α u (1 − u)
```

`u(x,y,t)` ∈ [0,1] is the fraction of people at a location who know the rumour. `D` is how far
information travels, `α` is how convincing it is. The logistic term saturates at `u = 1` — once
everyone at a point knows, it can't spread further there, only outward.

**Why this equation and not the heat equation.** Pure diffusion would let the rumour spread out and fade
— the peak drops as it widens, and eventually nothing is above the noise. The reaction term `αu(1−u)`
regenerates awareness locally, so instead of fading you get a travelling wave: a front that keeps its
shape and moves at constant speed until the whole population is converted. That's the behaviour models
of [social contagion](https://en.wikipedia.org/wiki/Social_contagion) and
[diffusion of innovations](https://en.wikipedia.org/wiki/Diffusion_of_innovations) are after, and it's
what makes Fisher-KPP the right choice here.

## Numerics

| | |
|---|---|
| Domain | Unit square [0,1]², Dirichlet `u = 0` on the boundary |
| Grid | 50 × 50 by default, `h = 1/49 ≈ 0.0204`; override with `-da_grid_x` / `-da_grid_y` |
| Spatial discretisation | 5-point Laplacian, 2nd-order central differences (`DMDA_STENCIL_STAR`, width 1) |
| Time integration | Explicit Runge–Kutta (`TSRK`), adaptive step |
| Decomposition | `DMDA` with `PETSC_DECIDE` in both directions; halo exchange via `DMGlobalToLocal` |
| Precision | `PetscScalar` = double |
| Initial condition | `u = 1` at the centre cell, `u = 0` everywhere else |
| Output | `.vts` every 5 steps into `files/` |

**Time step is limited by diffusion, not by the reaction.** The scheme is explicit, so it needs
`Δt ≲ h²/(4D)`. That bound moves a lot across the parameter range:

| Case | α | D | Wave speed `2√(αD)` | Δt limit | Steps for t = 10 |
|---|---:|---:|---:|---:|---:|
| Default | 0.5 | 0.001 | 0.045 | 0.104 | ~96 |
| Small-town gossip | 1.5 | 0.0001 | 0.024 | 1.04 | ~15 |
| Viral meme | 0.2 | 0.05 | 0.200 | 0.0021 | ~4800 |

The viral case needs 300× more timesteps than the gossip case for the same simulated time, entirely
because `D` is 500× larger. This is the standard trade-off with explicit time stepping on a diffusion
term, and the reason implicit methods exist. PETSc makes swapping one in a command-line flag — see
below.

## Checking it against theory

Fisher-KPP has a known travelling-wave solution moving at `c = 2√(αD)`. That gives a free correctness
check: measure how far the front has moved in the output and compare.

For the viral case, `c = 0.2`, so in 10 time units the front travels 2.0 — four times the 0.5
half-domain, meaning the wave crosses and saturates well before the run ends. For small-town gossip,
`c = 0.024` covers 0.24 in the same time, so it's still mid-domain when the run stops. Both match what
the `.vts` output shows.

## Build

Needs a configured PETSc install.

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=your-arch
make
```

The [`Makefile`](Makefile) pulls in PETSc's own `variables` and `rules`, so the include and link flags
match however PETSc was configured rather than being hardcoded.

## Run

**Small-town gossip** — sticky rumour, slow-moving:

```bash
mpiexec -n 4 ./rumour -alpha 1.5 -D 0.0001 -ts_max_time 10.0 -ts_monitor
```

**Viral meme** — weak persuasion, spreads fast:

```bash
mpiexec -n 4 ./rumour -alpha 0.2 -D 0.05 -ts_max_time 10.0 -ts_monitor
```

Everything PETSc exposes is available on the command line without recompiling:

```bash
-da_grid_x 200 -da_grid_y 200     # finer grid
-ts_type beuler                    # implicit, sidesteps the Δt limit above
-ts_adapt_type none -ts_dt 0.01    # fixed step
-log_view                          # PETSc performance summary
```

## Results

Output is a `.vts` series in `files/`. Open the directory as a database in ParaView or VisIt to animate
the front.

**Step 0** — one point of awareness:

![Frame 0](results/frame0.png)

**Step 115** — travelling wave has spread across the domain:

![Frame 115](results/frame115.png)

## Reading

- [Fisher's equation](https://en.wikipedia.org/wiki/Fisher%27s_equation) — the 1-D original and its
  travelling-wave solution
- [Bass diffusion model](https://en.wikipedia.org/wiki/Bass_diffusion_model) — the same idea in
  marketing, without the spatial term
- [Social contagion](https://en.wikipedia.org/wiki/Social_contagion) ·
  [Diffusion of innovations](https://en.wikipedia.org/wiki/Diffusion_of_innovations)

## Where this came from

Written March 2026 for *P2.4 — PETSc*, Master in High Performance Computing (ICTP / SISSA, Trieste).
It has been public since June 2026 in
[`prabhkodes/petsc`](https://github.com/prabhkodes/petsc), which also holds the introductory exercises
that led up to it — parallel vectors, sparse matrix assembly, KSP linear solves, a 2-D Poisson problem
on a `DMDA`, SNES nonlinear solves, and a pendulum ODE with `TS`. The course repo itself belongs to
SISSA and is private.

Two tidy-ups when moving the code here: an unused `MatStencil`/`PetscScalar` pair was removed from
`main`, and a stray `#include <sys/stat.h>` sitting in the middle of the file moved to the top. The
numerics are untouched.

## Layout

```
src/rumour.c    the solver — RHS function, VTK monitor, main
Makefile        uses PETSc's own build variables
results/        rendered frames
```
