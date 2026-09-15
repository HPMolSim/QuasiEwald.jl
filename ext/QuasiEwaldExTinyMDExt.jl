module QuasiEwaldExTinyMDExt

# Bridge between ExTinyMD's MD loop and QuasiEwald's framework-free plans.
#
# ## Why the interaction types live here, not in src/
#
# `MDSys`'s constructor requires `interactions::Vector{T_INTERACTION}` with
# `T_INTERACTION <: Tuple{ExTinyMD.AbstractInteraction, ExTinyMD.AbstractNeighborFinder}`.
# A struct's supertype is fixed where the struct is defined -- Julia has no
# mechanism for an extension, loaded later and conditionally, to retroactively
# add a supertype to an already-compiled type. `QuasiEwald`'s src/ does not
# depend on ExTinyMD at all (a weak dependency only), so nothing defined there
# can ever be a subtype of `ExTinyMD.AbstractInteraction` -- not "unless you
# remember an annotation", but structurally, in every build of the package.
# This is the same conclusion ParticleMeshEwald's extension documents, taken
# one step further: QuasiEwald has actual forces and drives `simulate!`
# (PME never could), so it needs real wrapper *types*, not just a method.
#
# `QuasiEwaldShortPlan`/`QuasiEwaldLongPlan`/`ZSorter` (src/types.jl) are the
# framework-free core: pure parameters (plus, for the long plan, no mass or
# scratch at all -- see its docstring). The three structs below are thin
# `ExTinyMD.AbstractInteraction`/`AbstractNeighborFinder` wrappers around one
# of those, holding only the MD-side scratch (gathered positions/charges, a
# force buffer) that a plan has no business owning. `QuasiEwald.jl`'s main
# module declares `QuasiEwaldShortInteraction`/`QuasiEwaldLongInteraction`/
# `SortingFinder` as dispatcher functions that call through to the
# constructors below via `Base.get_extension` once this extension has
# loaded, so `using QuasiEwald, ExTinyMD; QuasiEwaldShortInteraction(...)`
# keeps working exactly as it did before this package was decoupled.

using QuasiEwald, ExTinyMD, StaticArrays

# ----------------------------------------------------------------------------
# Gather helpers. Index convention (matching ExTinyMD's own adapter,
# ../ExTinyMD.jl/src/interactions/electrostatics/adapter.jl): `info.particle_info`
# is indexed by storage *slot*, `sys.atoms` by particle *id*. Positions are
# read in slot order and charges gathered to match, so plan-layer index `i`
# consistently means "slot i" and forces come back in slot order.
# ----------------------------------------------------------------------------

"Gather positions in slot order into `buf`, as SVector{3,T} (the plan's canonical AoS element)."
function gather_positions!(buf::Vector{SVector{3,T}}, info::ExTinyMD.SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        p = info.particle_info[i].position
        buf[i] = SVector{3,T}(p[1], p[2], p[3])
    end
    return buf
end

"Gather charges in slot order into `buf`, honouring the id/slot indirection."
function gather_charges!(buf::Vector{T}, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        buf[i] = sys.atoms[info.particle_info[i].id].charge
    end
    return buf
end

# A NoNeighborFinder carries no usable list, so the short-range plan falls
# back to its own O(n^2) pair loop in that case (see QuasiEwald.energy's
# docstring for QuasiEwaldShortPlan).
_finder_list(::ExTinyMD.NoNeighborFinder) = nothing
_finder_list(f) = f.neighbor_list

# ----------------------------------------------------------------------------
# Short-range wrapper
# ----------------------------------------------------------------------------

"""
    QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t)

`ExTinyMD.AbstractInteraction` wrapper around a [`QuasiEwald.QuasiEwaldShortPlan`](@ref).
Construct exactly as the pre-decoupling `QuasiEwaldShortInteraction` was
constructed; place `(interaction, finder)` in `sys.interactions` as before
(`finder` is an `ExTinyMD.CellListQ2D`, `CellListDirQ2D`, or `NoNeighborFinder`).
"""
struct QuasiEwaldShortInteraction{P,T} <: ExTinyMD.AbstractInteraction
    plan::P
    pos_scratch::Vector{SVector{3,T}}
    charge_scratch::Vector{T}
    force_buffer::Vector{SVector{3,T}}
end

function QuasiEwaldShortInteraction(plan::QuasiEwald.QuasiEwaldShortPlan{T}) where {T}
    n = plan.n_atoms
    return QuasiEwaldShortInteraction{typeof(plan),T}(
        plan, Vector{SVector{3,T}}(undef, n), Vector{T}(undef, n), Vector{SVector{3,T}}(undef, n))
end

QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t) =
    QuasiEwaldShortInteraction(QuasiEwald.QuasiEwaldShortPlan(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t))

function ExTinyMD.energy(inter::QuasiEwaldShortInteraction, neighborfinder,
                         sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    return QuasiEwald.energy(inter.plan, poses, charges; neighbor_list = _finder_list(neighborfinder))
end

function ExTinyMD.update_acceleration!(inter::QuasiEwaldShortInteraction, neighborfinder,
                                       sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    F = QuasiEwald.force!(inter.force_buffer, inter.plan, poses, charges;
                         neighbor_list = _finder_list(neighborfinder))
    @inbounds for i in eachindex(info.particle_info)
        m = sys.atoms[info.particle_info[i].id].mass
        f = F[i]
        info.particle_info[i].acceleration += ExTinyMD.Point(f[1] / m, f[2] / m, f[3] / m)
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Long-range wrapper
# ----------------------------------------------------------------------------

"""
    QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, k_c, rbe_p; Δk = ...)

`ExTinyMD.AbstractInteraction` wrapper around a [`QuasiEwald.QuasiEwaldLongPlan`](@ref).
Construct exactly as the pre-decoupling `QuasiEwaldLongInteraction` was
constructed; place `(interaction, finder)` in `sys.interactions` as before
(`finder` is a [`SortingFinder`](@ref)).
"""
struct QuasiEwaldLongInteraction{P,T} <: ExTinyMD.AbstractInteraction
    plan::P
    pos_scratch::Vector{SVector{3,T}}
    charge_scratch::Vector{T}
    force_buffer::Vector{SVector{3,T}}
end

function QuasiEwaldLongInteraction(plan::QuasiEwald.QuasiEwaldLongPlan{T}) where {T}
    n = plan.n_atoms
    return QuasiEwaldLongInteraction{typeof(plan),T}(
        plan, Vector{SVector{3,T}}(undef, n), Vector{T}(undef, n), Vector{SVector{3,T}}(undef, n))
end

QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, k_c, rbe_p; kwargs...) =
    QuasiEwaldLongInteraction(QuasiEwald.QuasiEwaldLongPlan(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, k_c, rbe_p; kwargs...))

# ----------------------------------------------------------------------------
# SortingFinder: an AbstractNeighborFinder wrapper around QuasiEwald.ZSorter,
# the plan-side half of the old SortingFinder (see src/types.jl).
# ----------------------------------------------------------------------------

"""
    SortingFinder(info::SimulationInfo)

`ExTinyMD.AbstractNeighborFinder` wrapper around a [`QuasiEwald.ZSorter`](@ref),
refreshed once per `simulate!` step (like every other finder) via
`ExTinyMD.update_finder!`. Pairs with [`QuasiEwaldLongInteraction`](@ref) in
`sys.interactions`.
"""
mutable struct SortingFinder{S} <: ExTinyMD.AbstractNeighborFinder
    sorter::S
end

function SortingFinder(info::ExTinyMD.SimulationInfo{T}) where {T}
    z_coords = [p_info.position[3] for p_info in info.particle_info]
    return SortingFinder(QuasiEwald.ZSorter(z_coords, sortperm(z_coords)))
end

function ExTinyMD.update_finder!(finder::SortingFinder, info::ExTinyMD.SimulationInfo)
    sorter = finder.sorter
    @inbounds for i in eachindex(info.particle_info)
        sorter.z_coords[i] = info.particle_info[i].position[3]
    end
    sortperm!(sorter.z_list, sorter.z_coords)
    return nothing
end

function ExTinyMD.energy(inter::QuasiEwaldLongInteraction, neighborfinder::SortingFinder,
                         sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    return QuasiEwald.energy(inter.plan, poses, charges; z_list = neighborfinder.sorter.z_list)
end

function ExTinyMD.update_acceleration!(inter::QuasiEwaldLongInteraction, neighborfinder::SortingFinder,
                                       sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    F = QuasiEwald.force!(inter.force_buffer, inter.plan, poses, charges;
                         z_list = neighborfinder.sorter.z_list)
    @inbounds for i in eachindex(info.particle_info)
        m = sys.atoms[info.particle_info[i].id].mass
        f = F[i]
        info.particle_info[i].acceleration += ExTinyMD.Point(f[1] / m, f[2] / m, f[3] / m)
    end
    return nothing
end

end
