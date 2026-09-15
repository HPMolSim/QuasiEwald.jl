module QuasiEwald

# these are packages to be used in this package
using LinearAlgebra, CellListMap, SpecialFunctions, GaussQuadrature, Distributions, Random, StaticArrays, StatsBase, Distributed
# ExTinyMD is imported narrowly (not a blanket `using`) so that this module's
# OWN `energy`/`force`/`force!` (the framework-free plan API, Task 3 of the
# decoupling phase) are genuinely new generic functions, not extensions of
# ExTinyMD's same-named `energy`. `using ExTinyMD` would bring `energy` into
# this module's unqualified scope bound to ExTinyMD's function, and Julia
# then requires (and this module does not want) every subsequent
# `function energy(...)` here to be a method of *that* function via
# `function ExTinyMD.energy(...)`. The names below are exactly the ones this
# module's still-ExTinyMD-coupled adapter code (to be moved to
# `ext/QuasiEwaldExTinyMDExt.jl` in the next task) uses unqualified; every
# use of `energy`/`update_acceleration!` on the ExTinyMD side is already
# spelled out fully as `ExTinyMD.energy`/`ExTinyMD.update_acceleration!`.
import ExTinyMD
using ExTinyMD: SimulationInfo, MDSys, CellListQ2D, update_finder!, Point

export IcmSys, GaussParameter, GreensElement, QuasiEwaldShortInteraction, QuasiEwaldLongInteraction, SortingFinder
export RBE_α, QuasiEwaldRbeInit
# Framework-free plans (Task 3 of the decoupling phase). `energy`/`force`/
# `force!` themselves are deliberately NOT exported -- see their docstrings.
export QuasiEwaldShortPlan, QuasiEwaldLongPlan, ZSorter, update_sorter!
export Gamma_1, Gamma_2, dz_Gamma_1, dz_Gamma_2, dz_Gamma_self_1, dz_Gamma_self_2
export IcmSysInit, IcmEnergy, IcmForce
export Gauss_int, Gauss_int_Tuple
export rbe_sampling
export Fsr_gauss_core, Fsr_point_core, Fsz_gauss_core, Fsz_point_core, Fsz_self_gauss_core, Fsz_self_point_core, QuasiEwald_Fs!, QuasiEwald_Fs_pair, QuasiEwald_Fs_self
export QuasiEwald_Fl!, force_long_total!, force_long_sampling!, force_direct_sum_k, force_long_k!, force_direct_sum_k0, force_k_sum_0, force_direct_sum_total
export energy_sum_total, energy_sum_sampling, QuasiEwald_El, Container, update_container!, direct_sum_total
export QuasiEwald_Es, QuaisEwald_Es_pair, QuaisEwald_Es_self, Es_gauss_core, Es_point_core
export RingAngles, nearest_angle_indice


include("types.jl")
include("init.jl")

include("tools/greens_functions.jl")
include("tools/Icm.jl")
include("tools/Gaussian_integrator.jl")
include("tools/Importance_sampling.jl")

include("energy/energy.jl")
include("energy/energy_long.jl")
include("energy/energy_short.jl")

include("force/force.jl")
include("force/force_long.jl")
include("force/force_short.jl") 

end


