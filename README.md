# QuasiEwald.jl

[![Build Status](https://github.com/ArrogantGao/QuasiEwald.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/QuasiEwald.jl/actions/workflows/CI.yml?query=branch%3Amain)

`QuasiEwald.jl` is a package written in `Julia`. 
It is an implementation of the algorithm Quasi Ewald Method, which used to calculate the electrostatic interaction in dielectric confined Quasi-2D charged systems for MD simulations, which has a linear complexity to the number of particles.

## Getting Started

`QuasiEwald.jl`'s core has **no dependency on ExTinyMD**: it takes plain
array-of-structs positions, a charge vector, and its own parameters, and
returns an energy or a force. `ExTinyMD` is only a *weak* dependency, used
solely to drive `simulate!`; it is loaded automatically (no configuration
needed) whenever you also `using ExTinyMD` in the same session, and its
absence in `[deps]` is why this package no longer requires it just to
compute an energy.

```julia
pkg> add https://github.com/ArrogantGao/QuasiEwald.jl
```

is enough for standalone use. Add `ExTinyMD` too if you want to drive an MD
loop (see [MD usage via ExTinyMD](#md-usage-via-extinymd) below).

> **`pkg> add QuasiEwald` will not get you this version.** The General
> registry's newest QuasiEwald is 0.2.1, from before the decoupling, and this
> package cannot be registered while its `Project.toml` carries the
> `[sources]` override described [further down](#a-note-for-anyone-pinning-this-package-against-a-sibling-extinymd-checkout).
> Install from the repository URL (above), or `pkg> dev` a local checkout,
> until that override can be removed.

### Standalone usage (no ExTinyMD)

The core API is a **plan** (pure parameters plus solver scratch, no
ExTinyMD type anywhere) queried against plain arrays:

```julia
using QuasiEwald, StaticArrays

n = 100
L = (100.0, 100.0, 10.0)     # (Lx, Ly, Lz); z is the non-periodic, confined axis
poses = [SVector(L[1]*rand(), L[2]*rand(), 1.0 + 8.0*rand()) for _ in 1:n]
charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

γ_1, γ_2 = 0.4, -0.5          # dielectric mismatch at the z = 0 / z = Lz boundaries
ϵ_0 = 1.0
accuracy = 1e-4
α = 1.0
r_c = 4.5                    # r_c must be < min(Lx, Ly) / 2 = 50.0 here
k_c = sqrt(-4 * α * log(accuracy))
n_t = 100                    # 30 suffices for energies but not for the
                             # short-range z-force -- see the note below

short_plan = QuasiEwaldShortPlan(γ_1, γ_2, ϵ_0, L, false, accuracy, α, n, r_c, n_t)
long_plan  = QuasiEwaldLongPlan(γ_1, γ_2, ϵ_0, L, false, accuracy, α, n, k_c, 0)

E = QuasiEwald.energy(short_plan, poses, charges) + QuasiEwald.energy(long_plan, poses, charges)
F = QuasiEwald.force(short_plan, poses, charges) .+ QuasiEwald.force(long_plan, poses, charges)
```

A few things worth knowing:

- **`poses` may be `Vector{SVector{3,T}}`, `Vector{NTuple{3,T}}`, or ExTinyMD's `Vector{Point{3,T}}`** --
  the kernels only ever index `p[1]`/`p[2]`/`p[3]`, never rely on
  vector-vector arithmetic on the elements you pass in, so no conversion
  layer is needed either way. `SVector{3,T}` is the canonical/recommended
  choice.
- **`QuasiEwald.energy`, `QuasiEwald.force` and `QuasiEwald.force!` are
  defined but deliberately not exported** -- always call them qualified.
  (If every electrostatics package in this family exported its own
  `energy`, `using QuasiEwald, SoEwald2D` would make the bare name
  ambiguous.) The plan **constructors** (`QuasiEwaldShortPlan`,
  `QuasiEwaldLongPlan`, `ZSorter`, and the ICM helpers `IcmSys`/
  `IcmSysInit`/`IcmEnergy`/`IcmForce`) are exported as usual.
- **`r_c` must be strictly less than `min(Lx, Ly) / 2`.** The geometry is
  quasi-2D: `x`/`y` are periodic, `z` is not (it is the confined,
  dielectric-bounded axis). The short-range cutoff only ever needs to see
  at most one periodic image per axis; anything ≥ half the box breaks that
  and is not supported. `QuasiEwaldShortPlan` throws an `ArgumentError`
  naming the offending `r_c`/`Lx`/`Ly` rather than returning a silently
  multiply-counted sum.
- Neither `poses` nor `charges` is ever mutated by a query.
- **`accuracy` and `n_t` must be tightened together when you care about
  forces.** `accuracy` sets where the k-space integrals are truncated; `n_t`
  sets the Gauss quadrature order over that interval. Raising `accuracy`
  alone makes the short-range **z**-force *worse*, not better, because a
  tighter truncation widens the interval that a fixed-order rule has to
  cover. Measured max `|F_z - (-dE/dz)|` for a 6-particle configuration:

  | `n_t` \ `accuracy` | `1e-4` | `1e-6` | `1e-8` | `1e-10` |
  |---|---|---|---|---|
  | 30  | 5.7e-6 | 7.5e-6 | 2.0e-5 | 3.8e-5 |
  | 60  | 6.3e-6 | 6.5e-8 | 4.1e-9 | 1.5e-8 |
  | 100 | 6.3e-6 | 6.5e-8 | 8.6e-10 | 1.4e-10 |
  | 400 | 6.3e-6 | 6.5e-8 | 8.6e-10 | 1.4e-10 |

  With `n_t` converged (≳ 100 here) the force error tracks `accuracy` at
  roughly `60 × accuracy`. `n_t = 30` is enough for energies but not for the
  z-derivative. The energy and the in-plane force components are far less
  sensitive to `n_t` than the z-force is.
- `QuasiEwald.energy(::QuasiEwaldShortPlan, poses, charges; ...)` (and
  likewise `force`/`force!`) accepts an optional `neighbor_list =` keyword
  (candidate `(i, j, ...)` pairs -- e.g. a `CellListMap` neighbor list you
  already maintain); the true in-plane distance is always recomputed from
  `poses`, so a supplied list's own reported distance is ignored. With none
  given, every pair is tested directly (`O(n_atoms^2)`) -- this plan does
  not own a persistent cell list of its own. It also accepts
  `single_mode = false`, which when set to `true` drops the
  Gaussian-screened part of each pair and self term; the default is the full
  short-range sum, and it is the only setting that pairs with
  `QuasiEwaldLongPlan` to give the correct total.
  `QuasiEwald.energy(::QuasiEwaldLongPlan, poses, charges; ...)` similarly
  accepts `z_list =` (a z-sort you already have, e.g. from a `ZSorter` you
  keep across calls -- `ZSorter(poses)` or `ZSorter(z_coords)`, refreshed in
  place with `update_sorter!`) and otherwise sorts fresh with `sortperm` each
  call.

### MD usage via ExTinyMD

Placing quasi-2D electrostatics in `sys.interactions` (so `simulate!`
drives it) needs `ExTinyMD.AbstractInteraction`/`AbstractNeighborFinder`
wrapper types, which can only be *defined* once ExTinyMD exists -- they
live in this package's `ExTinyMD` extension (`ext/QuasiEwaldExTinyMDExt.jl`),
loaded automatically the moment both packages are `using`'d together.
`QuasiEwaldShortInteraction`, `QuasiEwaldLongInteraction` and
`SortingFinder` are constructed exactly as before this package was
decoupled from ExTinyMD -- the wrapper just holds a
`QuasiEwaldShortPlan`/`QuasiEwaldLongPlan`/`ZSorter` internally now:

```julia
using ExTinyMD, QuasiEwald

intershort = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, r_c, n_t)
short_finder = CellListQ2D(info, r_c + 1.0, boundary, 100)
interlong = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, k_c, rbe_p)
long_finder = SortingFinder(info)

sys = MDSys(
    n_atoms = n_atoms, atoms = atoms, boundary = boundary,
    interactions = [(intershort, short_finder), (interlong, long_finder)],
    loggers = loggers, simulator = simulator,
)
simulate!(simulator, sys, info, n_steps)
```

Calling `QuasiEwaldShortInteraction`/`QuasiEwaldLongInteraction`/
`SortingFinder` before `using ExTinyMD` raises a clear error naming the
extension, rather than an `UndefVarError`. If ExTinyMD *is* loaded and the
error still appears, it says so and points at `Base.retry_load_extensions()`
-- that case means the extension itself failed to precompile, not that
`using ExTinyMD` is missing.

`short_finder` must be an `ExTinyMD.CellListQ2D` or `CellListDirQ2D` (both
select candidate pairs by *in-plane* distance), or a `NoNeighborFinder` to
fall back to the plan's own `O(n_atoms^2)` pair loop. A 3-D finder
(`CellList3D`, `CellListDir3D`) is refused with an `ArgumentError`: its list
is built from 3-D distances and so omits pairs that are close in plane but
far apart in `z`, which would make the quasi-2D short-range sum silently too
small.

#### Breaking change: the three exported names are functions, not types

Before the decoupling, `QuasiEwaldShortInteraction`,
`QuasiEwaldLongInteraction` and `SortingFinder` were `struct`s defined in
`src/`. They cannot be, any more: a struct's supertype is fixed where the
struct is defined, and `src/` has no ExTinyMD dependency to name
`ExTinyMD.AbstractInteraction` with. The real types now live in the
extension, and the three exported names are *dispatcher functions* that
forward to them. Consequently:

- **Construction works unchanged.** Every constructor call above behaves
  exactly as it did before, and the result still
  `isa ExTinyMD.AbstractInteraction`.
- **Type-position uses do not work.** `x isa QuasiEwaldShortInteraction`, an
  `::QuasiEwaldShortInteraction` argument annotation, a
  `Vector{QuasiEwaldShortInteraction}` element type, and dispatching your own
  method on one all now raise
  `TypeError: in isa, expected Type, got a value of type typeof(QuasiEwaldShortInteraction)`.

If you need the type itself, reach into the extension module:

```julia
using QuasiEwald, ExTinyMD
ext = Base.get_extension(QuasiEwald, :QuasiEwaldExTinyMDExt)

x isa ext.QuasiEwaldShortInteraction                  # works
f(x::ext.QuasiEwaldLongInteraction) = ...             # works
```

The type cannot also be re-exported from `QuasiEwald` under the same name,
because that exported name is precisely what lets the constructor call
resolve without ExTinyMD being a hard dependency.

### A note for anyone pinning this package against a sibling ExTinyMD checkout

This package's `Project.toml` currently carries a `[sources]` override
pinning `ExTinyMD` to its GitHub `main` branch, because ExTinyMD 0.3 (the
version this package requires) is not yet on the General registry (only
0.2.7 is). That override is temporary and tracked centrally across every
downstream package in
[`ExTinyMD.jl`'s downstream-decoupling design doc, §6b](https://github.com/HPMolSim/ExTinyMD.jl/blob/main/docs/superpowers/specs/2026-09-15-downstream-decoupling-design.md);
General's automerge refuses any package whose `Project.toml` contains
`[sources]`, so this package cannot be tagged/registered until ExTinyMD 0.3
is registered and this override is removed.

Here is an simple example, which will calculate the interaction between two paricle confined by dielectric substrate of different dielectric permittivity.
```julia
using Plots, ExTinyMD, QuasiEwald

begin
    n_atoms = 2
    L = 180.0
    boundary = ExTinyMD.Q2dBoundary(L, L, 10.0) 

    atoms = [Atom(type = 1, mass = 1.0, charge = 1.0), Atom(type = 2, mass = 1.0, charge = -1.0)]

    sys = MDSys(
                n_atoms = n_atoms,
                atoms = atoms,
                boundary = boundary,
                interactions = [(NoInteraction(), NoNeighborFinder(n_atoms))],
                loggers = [TrajectionLogger(step = 100, output = false)],
                simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
            )

    Force_x = Vector{Vector{Float64}}()
    coord_1 = Point(50.0, 50.0, 1.0)
    X = 0.1:0.1:40.0

    for (γ_1, γ_2) in [(0.0, 0.0), (0.95, 0.95), (-0.95, -0.95), (10.0, 10.0), (-10.0, -10.0)]
        ϵ_0 = 1.0
        n_t = 30

        accuracy = 1e-4
        α = 1.0
        r_c = 15.0
        k_c = sqrt(-4 * α * log(accuracy))

        force_x = Vector{Float64}()

        for x_2 in X
            info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, 10.0), boundary; min_r = 1.0, temp = 1.0)
            coord_2 = Point(50.0 + x_2, 50.0, 1.01)

            info.particle_info[1].position = coord_1
            info.particle_info[2].position = coord_2

            sortz = SortingFinder(info)
            cellq2d = CellListQ2D(info, r_c, boundary, 1)
            interaction_short = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, r_c, n_t)
            interaction_long = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, k_c, 0)
    
            info.particle_info[1].acceleration = Point(0.0, 0.0, 0.0)
            info.particle_info[2].acceleration = Point(0.0, 0.0, 0.0)
            ExTinyMD.update_acceleration!(interaction_short, cellq2d, sys, info)
            ExTinyMD.update_acceleration!(interaction_long, sortz, sys, info)
            push!(force_x, info.particle_info[1].acceleration[1])
        end
        push!(Force_x, force_x)
    end

    plot(dpi = 300, size = (800, 600), legend = :topright, xlabel = "x", ylabel = "force_x")
    for (γ, force_x) in zip([0.0, 0.95, -0.95, 10.0, -10.0], Force_x)
        plot!(X, force_x, label = "γ = " * string(γ), ylim = [-0.06, 0.06])
    end
    savefig("force_x.png")
end
```
To run this script, you can simple type:
```
julia ./example/force/force_pair.jl
```
The result is shown below:

![Force in x direction](./examples/force/force_x.png)


Here is another example, which shows how to simulate a dielectric confined charged system via `ExTinyMD.jl` and `QuasiEwald.jl`.
```julia
using ExTinyMD, QuasiEwald

begin
    n_atoms = 436
    n_atoms = Int64(round(n_atoms))
    L_x = 100.0
    L_y = 100.0
    L_z = 50.0
    L = (L_x, L_y, L_z)
    boundary = Q2dBoundary(L_x, L_y, L_z)
    atoms = Vector{Atom{Float64}}()

    for i in 1:218
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in 219:436
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    (γ_1, γ_2) = (0.95, -0.95)

    info = SimulationInfo(n_atoms, atoms, (0.0, L_x, 0.0, L_y, 0.5, L_z - 0.5), boundary; min_r = 2.0, temp = 1.0)

    ϵ_0 = 1.0

    accuracy = 1e-4
    α = 1.0
    k_c = sqrt(- 4 * α * log(accuracy))
    r_c = (α * accuracy)^(-1/3) / 2
    n_t = 30
    rbe_p = 50

    intershort = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, r_c, n_t)
    short_finder = CellListQ2D(info, r_c + 1.0, boundary, 100)
    interlong = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, k_c, rbe_p)
    long_finder = SortingFinder(info)

    interactions = [
        (LennardJones(), CellList3D(info, 4.5, boundary, 100)),
        (SubLennardJones(0.0, L_z; cutoff = 0.5, σ = 0.5), SubNeighborFinder(1.0, info, 0.0, L_z)), 
        (intershort, short_finder),
        (interlong, long_finder)
        ]

    loggers = [TempartureLogger(100, output = true), TrajectionLogger(step = 100, output = true)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    simulate!(simulator, sys, info, 100000)
end
```
which will simulate 436 changed particle confined by substrates with $\gamma = 0.95$ for $10^6$ steps, and the trajection will be recorded.

## Questions and Contributions

Please open an [issue](https://github.com/ArrogantGao/QuasiEwald.jl/issues)
if you encounter any problems, or have any feature requests.

It is also welcomed for any suggestions about the issues marked as `enhancement`, please let us know if you have any idea about them.