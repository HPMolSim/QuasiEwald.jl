module QuasiEwald

# these are packages to be used in this package. ExTinyMD is NOT here: it is a
# [weakdeps] entry only (see Project.toml and ext/QuasiEwaldExTinyMDExt.jl).
# This module has no ExTinyMD dependency at all -- everything below is plain
# arrays and numbers.
# `CellListMap` and `Distributions` used to be listed here and were both
# dead: this package never calls CellListMap at all (the only mention left is
# a docstring, describing what a caller may pass as `neighbor_list`), and the
# one sampling call in tools/Importance_sampling.jl -- `sample(K_set,
# ProbabilityWeights(Prob), 1000)` -- comes from StatsBase, not Distributions.
# `Distributed` IS live (`@distributed (+)` in force/force_long.jl).
using LinearAlgebra, SpecialFunctions, GaussQuadrature, Random, StaticArrays, StatsBase, Distributed

export IcmSys, GaussParameter, GreensElement
export RBE_α
# QuasiEwaldShortInteraction / QuasiEwaldLongInteraction / SortingFinder are
# ExTinyMD.AbstractInteraction/AbstractNeighborFinder wrappers that can only
# be *defined* once ExTinyMD exists (a struct's supertype is fixed where the
# struct is defined; see ext/QuasiEwaldExTinyMDExt.jl's module docstring).
# Since this module cannot know ExTinyMD's types at all, these three names
# are declared here as plain dispatcher functions -- exported so `using
# QuasiEwald, ExTinyMD` keeps working exactly as before -- that hand off to
# the extension's real constructors via `Base.get_extension`. This is the
# standard pattern for a package extension that must introduce a brand new
# type (as opposed to adding a method to an existing function): a struct
# definition cannot be written as `struct QuasiEwald.Foo ... end` (dot-
# qualified struct definitions are not legal Julia), so the type itself has
# to live in the extension, and this stub is what makes `QuasiEwaldShortInteraction(...)`
# resolve to it without every caller needing to know that.
export QuasiEwaldShortInteraction, QuasiEwaldLongInteraction, SortingFinder
# Framework-free plans (Task 3 of the decoupling phase). `energy`/`force`/
# `force!` themselves are deliberately NOT exported -- see their docstrings.
export QuasiEwaldShortPlan, QuasiEwaldLongPlan, ZSorter, update_sorter!
export Gamma_1, Gamma_2, dz_Gamma_1, dz_Gamma_2, dz_Gamma_self_1, dz_Gamma_self_2
export IcmSysInit, IcmEnergy, IcmForce
export Gauss_int, Gauss_int_Tuple
export rbe_sampling
export Fsr_gauss_core, Fsr_point_core, Fsz_gauss_core, Fsz_point_core, Fsz_self_gauss_core, Fsz_self_point_core, QuasiEwald_Fs_pair, QuasiEwald_Fs_self
export force_long_total!, force_long_sampling!, force_direct_sum_k, force_long_k!, force_direct_sum_k0, force_k_sum_0, force_direct_sum_total
export energy_sum_total, energy_sum_sampling, Container, update_container!, direct_sum_total
export QuaisEwald_Es_pair, QuaisEwald_Es_self, Es_gauss_core, Es_point_core
export RingAngles, nearest_angle_indice

"""
    QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t)
    QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, k_c, rbe_p; Δk = ...)
    SortingFinder(info)

`ExTinyMD.AbstractInteraction`/`AbstractNeighborFinder` wrappers around a
[`QuasiEwaldShortPlan`](@ref)/[`QuasiEwaldLongPlan`](@ref)/[`ZSorter`](@ref),
for use in `sys.interactions`. Defined by `ext/QuasiEwaldExTinyMDExt.jl`,
loaded automatically once `using ExTinyMD` has also been done -- calling
these before that raises an informative error rather than a cryptic
`UndefVarError`. For standalone use (no ExTinyMD, no MDSys), construct a
plan directly and query it with `QuasiEwald.energy`/`force`/`force!`.

!!! warning "Breaking change: these three names are functions, not types"
    Before this package was decoupled from ExTinyMD, each of these was a
    `struct`. They are now *dispatcher functions* that forward to the real
    constructors in the extension, so:

    * **Construction works unchanged.** `QuasiEwaldShortInteraction(γ_1, γ_2,
      ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t)` returns exactly the
      wrapper it always did, and it still `isa ExTinyMD.AbstractInteraction`.
    * **Type-position uses do not work.** `x isa QuasiEwaldShortInteraction`,
      an `::QuasiEwaldShortInteraction` annotation, a
      `Vector{QuasiEwaldShortInteraction}` element type, and dispatching a
      method on one all now raise
      `TypeError: in isa, expected Type, got a value of type typeof(QuasiEwaldShortInteraction)`.

    If you need the type itself, reach into the extension module:

    ```julia
    using QuasiEwald, ExTinyMD
    ext = Base.get_extension(QuasiEwald, :QuasiEwaldExTinyMDExt)
    x isa ext.QuasiEwaldShortInteraction        # works
    ```

    The type cannot be re-exported from this module under the same name,
    because the exported name is what makes the constructor call resolve
    without ExTinyMD being a hard dependency.
"""
QuasiEwaldShortInteraction, QuasiEwaldLongInteraction, SortingFinder

# ExTinyMD's PkgId, for telling "ExTinyMD was never loaded" apart from
# "ExTinyMD is loaded but the extension did not come up". `Base.get_extension`
# returns `nothing` in both cases, and reporting the second as the first sends
# the user off to `using ExTinyMD` -- which they have already done -- while the
# actual precompilation error has scrolled off the screen.
const _EXTINYMD_PKGID = Base.PkgId(Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD")

for name in (:QuasiEwaldShortInteraction, :QuasiEwaldLongInteraction, :SortingFinder)
    @eval function $name(args...; kwargs...)
        ext = Base.get_extension(QuasiEwald, :QuasiEwaldExTinyMDExt)
        if ext === nothing
            if haskey(Base.loaded_modules, _EXTINYMD_PKGID)
                error($(string(name)) * ": ExTinyMD is loaded, but QuasiEwald's " *
                      "QuasiEwaldExTinyMDExt extension failed to load -- so `using ExTinyMD` " *
                      "is not what is missing. This is almost always a precompilation error " *
                      "inside the extension, reported as a `Error: Error during loading of " *
                      "extension QuasiEwaldExTinyMDExt of QuasiEwald` warning that has since " *
                      "scrolled past. Run `Base.retry_load_extensions()` to reproduce it, and " *
                      "check for a version conflict between the loaded ExTinyMD and this " *
                      "package's [compat] bound.")
            else
                error($(string(name)) * " requires ExTinyMD to be loaded (`using ExTinyMD`) -- " *
                      "it is defined by the QuasiEwaldExTinyMDExt package extension. For standalone " *
                      "use, construct a QuasiEwaldShortPlan/QuasiEwaldLongPlan directly instead.")
            end
        end
        return getfield(ext, $(QuoteNode(name)))(args...; kwargs...)
    end
end

include("types.jl")
include("init.jl")

include("tools/greens_functions.jl")
include("tools/Icm.jl")
include("tools/Gaussian_integrator.jl")
include("tools/Importance_sampling.jl")

include("energy/energy_long.jl")
include("energy/energy_short.jl")

include("force/force_long.jl")
include("force/force_short.jl")

end
