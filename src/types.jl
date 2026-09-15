# Nearest-image displacement, `dx - L*round(dx/L)`. This package used to rely on
# ExTinyMD's `position_checkQ2D`, which scans `m ∈ -1:1` and returns the *first*
# periodic image inside the cutoff (or a sentinel all-zero triple if none is
# found, forcing every caller to guard with `iszero`) -- correct only while
# `r_c < L/2`. `_wrap`/`_min_image_q2d` are the true minimum image for any cutoff,
# and callers test the returned `ρ_sq` explicitly instead of a sentinel.
@inline _wrap(dx::T, L::T) where {T} = dx - L * round(dx / L)

"""
    _min_image_q2d(pos_i, pos_j, L) -> (coord_i, coord_j, ρ_sq)

Quasi-2D nearest-image helper: x and y wrap under `L`, z is left as a plain
difference (it is not periodic in this geometry). Returns `pos_i` shifted to
its nearest in-plane image of `pos_j`, `pos_j` unchanged, and their squared
in-plane distance `ρ_sq`. `pos_i`/`pos_j` need only support `p[1]`, `p[2]`,
`p[3]` indexing (an `SVector{3,T}`, `NTuple{3,T}` or ExTinyMD's `Point{3,T}`
all qualify).
"""
@inline function _min_image_q2d(pos_i, pos_j, L::NTuple{3, T}) where {T}
    dx = _wrap(T(pos_i[1]) - T(pos_j[1]), L[1])
    dy = _wrap(T(pos_i[2]) - T(pos_j[2]), L[2])
    ρ_sq = dx^2 + dy^2
    coord_i = SVector{3, T}(T(pos_j[1]) + dx, T(pos_j[2]) + dy, T(pos_i[3]))
    coord_j = SVector{3, T}(T(pos_j[1]), T(pos_j[2]), T(pos_j[3]))
    return coord_i, coord_j, ρ_sq
end

struct IcmSys{T, R}
    γ::NTuple{2, T} # (γ_up, γ_down)
    L::NTuple{3, T} # (Lx, Ly, Lz)
    N_real::R
    N_img::R
end

IcmSys(γ::NTuple{2, Float64}, L::NTuple{3, Float64}, N_real::Int, N_img::Int) = IcmSys{Float64, Int}(γ, L, N_real, N_img)

struct GaussParameter{T}
    sw::Vector{NTuple{2, T}}
end

GaussParameter(Step::Int) = GaussParameter{Float64}([tuple(legendre(Step)[1][i], legendre(Step)[2][i]) for i in 1:Step])

struct GreensElement{T}
    γ_1::T
    γ_2::T
    ρ::T
    a::NTuple{4, T}
    b::NTuple{4, T}
    sign_a::NTuple{4, T}
    L_z::T
    α::T
    k_f1::NTuple{4, T}
    k_f2::NTuple{4, T}
end


GreensElement(γ_1::T, γ_2::T, L_z::T, α::T) where T = GreensElement{T}(γ_1, γ_2, zero(T), (zero(T), zero(T), zero(T), zero(T)), (zero(T), zero(T), zero(T), zero(T)), (zero(T), -one(T), one(T), zero(T)), L_z, α, (zero(T), zero(T), zero(T), zero(T)), (zero(T), zero(T), zero(T), zero(T)))

function GreensElement(γ_1::T, γ_2::T, z_i::T, L_z::T, α::T, accuracy::T) where T
    ρ = zero(T)
    z_n = zero(T)
    z_p = 2 * z_i

    a = (z_n, z_p, 2 * L_z - z_p, 2 * L_z - z_n)
    b = (1.0, γ_1, γ_2, γ_1 * γ_2)

    sign_a = (zero(T), -one(T), one(T), zero(T))
    k_f1 = sqrt.(4 * α^2 .* a.^2 .- 4 * α * log(accuracy)) .- 2 * α .* a
    k_f2 = - log(accuracy) ./ (2 * L_z .+ a)

    return GreensElement{T}(γ_1, γ_2, ρ, a, b, sign_a, L_z, α, k_f1, k_f2)
end

function GreensElement(γ_1::T, γ_2::T, z_i::T, z_j::T, ρ::T, L_z::T, α::T, accuracy::T) where T
    z_n = abs(z_i - z_j)
    z_p = z_i + z_j

    a = (z_n, z_p, 2 * L_z - z_p, 2 * L_z - z_n)
    b = (1.0, γ_1, γ_2, γ_1 * γ_2)

    sign_a = (-sign(z_i - z_j), -one(T), one(T), sign(z_i - z_j))
    k_f1 = sqrt.(4 * α^2 .* a.^2 .- 4 * α * log(accuracy)) .- 2 * α .* a
    k_f2 = - log(accuracy) ./ (2 * L_z .+ a)

    return GreensElement{T}(γ_1, γ_2, ρ, a, b, sign_a, L_z, α, k_f1, k_f2)
end

struct RingAngles{T}
    ring_angles::Vector{T}
    sectors_sum::Vector{T}
end

function RingAngles(k_0::T, L_x::T, L_y::T, L_z::T, α::T, k_c::T, Δk::T) where{T}
    ring_angles = Vector{T}()
    nx_max = ceil(Int, k_0 * L_x / 2π) + 2
    ny_max = ceil(Int, k_0 * L_y / 2π) + 2

    push!(ring_angles, -π)
    for nx in - nx_max : nx_max
        for ny in - ny_max : ny_max
            k_x = nx * 2π / L_x
            k_y = ny * 2π / L_y
            k = sqrt(k_x^2 + k_y^2)
            if abs(k - k_0) < Δk
                θ = atan(k_y, k_x)
                push!(ring_angles, θ)
            end
        end
    end
    sort!(ring_angles)

    sectors_sum = zeros(T, length(ring_angles) - 1)

    nxc_max = ceil(Int, k_c * L_x / 2π) + 1
    nyc_max = ceil(Int, k_c * L_y / 2π) + 1

    for nx in - nxc_max : nxc_max
        for ny in - nyc_max : nyc_max
            k_x = nx * 2π / L_x
            k_y = ny * 2π / L_y
            k = sqrt(k_x^2 + k_y^2)
            if 0 < k <= k_c 
                nearest_id = nearest_angle_indice(k_x, k_y, ring_angles)
                sectors_sum[nearest_id] += exp(- k^2 / (4 * α)) / (exp(- 2 * (k - k_0) * L_z) - 1)
            end
        end
    end

    return RingAngles{T}(ring_angles, sectors_sum)
end

function RingAngles(k_0::T) where{T}
    return RingAngles{T}(Vector{T}(), Vector{T}())
end

function nearest_angle_indice(k_x::T, k_y::T, ring_angles::Vector{T}) where{T}
    θ_k = atan(k_y, k_x)
    θ_ring = ring_angles

    low = 1
    high = length(θ_ring)
    closest = θ_ring[low]
    id = 1

    while low <= high
        mid = low + (high - low) ÷ 2
        if θ_ring[mid] == θ_k
            id = mid
            break
        elseif θ_ring[mid] < θ_k
            low = mid + 1
        else
            high = mid - 1
        end

        if abs(θ_ring[mid] - θ_k) < abs(closest - θ_k)
            closest = θ_ring[mid]
            id = mid
        end
    end

    if id == length(θ_ring)
        id = 1
    end

    return id
end

# ============================================================================
# Framework-free plans.
#
# `QuasiEwaldShortInteraction`/`QuasiEwaldLongInteraction`/`SortingFinder`
# (the dispatcher stubs declared in QuasiEwald.jl) used to be defined here as
# structs wearing an `ExTinyMD.AbstractInteraction`/`AbstractNeighborFinder`
# supertype directly. That supertype is fixed at struct definition and
# cannot be retrofitted by an extension (see the phase spec's §4.3a), which
# is fundamentally incompatible with this module having no ExTinyMD
# dependency at all -- so those struct definitions, and everything MDSys/
# SimulationInfo-shaped that went with them, now live in
# `ext/QuasiEwaldExTinyMDExt.jl`, defined only once ExTinyMD is loaded.
#
# `QuasiEwaldShortPlan`/`QuasiEwaldLongPlan`/`ZSorter` below are the
# framework-free replacement: the same physics, no ExTinyMD anywhere,
# constructed and queried from plain arrays via `QuasiEwald.energy`/`force`/
# `force!` (defined alongside the rest of the energy/force machinery in
# energy/energy_short.jl, energy/energy_long.jl, force/force_short.jl,
# force/force_long.jl). `mass`/`acceleration` do not appear on the long
# plan: a framework-free solver returns forces and leaves mass-division to
# the caller, exactly as ExTinyMD's own electrostatics adapter does
# (../ExTinyMD.jl/src/interactions/electrostatics/adapter.jl) -- the
# extension's `update_acceleration!` does that division.
# ============================================================================

"""
    QuasiEwaldShortPlan(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t)

Framework-free short-range (real-space) plan for the quasi-2D Ewald method.
Pure parameters -- no ExTinyMD dependency, nothing MD-specific. Query with
[`QuasiEwald.energy`](@ref), [`QuasiEwald.force`](@ref) or
[`QuasiEwald.force!`](@ref) against plain array-of-structs positions
(`Vector{SVector{3,T}}` canonical, but anything supporting `p[1]`/`p[2]`/
`p[3]` indexing works) and a plain charge vector.

`r_c` must satisfy `r_c < min(Lx, Ly) / 2`; anything else throws an
`ArgumentError` (the short-range sum uses the single nearest in-plane
periodic image, which is only correct below half the box).
"""
struct QuasiEwaldShortPlan{T, TI}
    γ_1::T
    γ_2::T
    ϵ_0::T
    L::NTuple{3, T}
    rbe::Bool
    accuracy::T
    α::T
    n_atoms::TI

    r_c::T
    n_t::TI
    gauss_para::GaussParameter{T}
end

function QuasiEwaldShortPlan(γ_1::T, γ_2::T, ϵ_0::T, L::NTuple{3, T}, rbe::Bool, accuracy::T, α::T, n_atoms::TI, r_c::T, n_t::TI) where {T<:Number, TI<:Integer}
    # `_min_image_q2d` returns the single nearest in-plane image, which is the
    # only image inside the cutoff exactly when `r_c < min(Lx, Ly) / 2`. At or
    # beyond half the box a second image is also within `r_c` and the
    # short-range sum silently multiply-counts it: no exception, no warning,
    # just a wrong number. This is the only place that can catch it for a
    # standalone (no-ExTinyMD, no-CellListMap) caller, since nothing else in
    # this package ever looks at the unit cell.
    if !(r_c < min(L[1], L[2]) / 2)
        throw(ArgumentError(
            "QuasiEwaldShortPlan requires r_c < min(Lx, Ly) / 2, but got " *
            "r_c = $r_c with (Lx, Ly) = ($(L[1]), $(L[2])), i.e. " *
            "min(Lx, Ly) / 2 = $(min(L[1], L[2]) / 2). The quasi-2D " *
            "short-range sum uses the single nearest in-plane periodic image, " *
            "which is only the whole story below half the box; at or above it " *
            "the sum multiply-counts images and the result is silently wrong. " *
            "Reduce r_c, or enlarge Lx/Ly."))
    end
    return QuasiEwaldShortPlan{T, TI}(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t, GaussParameter(n_t))
end

"""
    QuasiEwaldLongPlan(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, k_c, rbe_p; Δk = ...)

Framework-free long-range (reciprocal-space) plan for the quasi-2D Ewald
method. Same parameters as the old `QuasiEwaldLongInteraction`, minus
`mass`/`coords`/`acceleration`: [`QuasiEwald.force!`](@ref) returns a force,
not an acceleration, so this plan carries no notion of mass at all. Query
with [`QuasiEwald.energy`](@ref)/[`QuasiEwald.force`](@ref)/
[`QuasiEwald.force!`](@ref); pass `z_list =` to reuse a z-sort already
computed (e.g. by [`ZSorter`](@ref)), or omit it to have the query compute
its own via `sortperm`.
"""
struct QuasiEwaldLongPlan{T, TI}
    γ_1::T
    γ_2::T
    ϵ_0::T
    L::NTuple{3, T}
    rbe::Bool
    accuracy::T
    α::T
    n_atoms::TI

    k_c::T
    rbe_p::TI
    sum_k::T
    K_set::Vector{NTuple{3, T}}

    k_0::T
    ringangles::RingAngles{T}
end

function QuasiEwaldLongPlan(γ_1::T, γ_2::T, ϵ_0::T, L::NTuple{3, T}, rbe::Bool, accuracy::T, α::T, n_atoms::TI, k_c::T, rbe_p::TI; Δk::T = π / sqrt(L[1] * L[2])) where {T<:Number, TI<:Integer}
    K_set, sum_k = rbe_sampling(L, α, accuracy)

    if γ_1 * γ_2 ≥ one(T)
        k_0 = log(γ_1 * γ_2) / (2 * L[3])
        ringangles = RingAngles(k_0, L[1], L[2], L[3], α, k_c, Δk)
    else
        k_0 = zero(T)
        ringangles = RingAngles(k_0)
    end

    return QuasiEwaldLongPlan{T, TI}(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, k_c, rbe_p, sum_k, K_set, k_0, ringangles)
end

"""
    ZSorter(poses) -> ZSorter
    ZSorter(z_coords::Vector{T}) -> ZSorter

Plan-side z-sorter: the framework-free half of the old `SortingFinder`
(itself an `ExTinyMD.AbstractNeighborFinder`, which cannot live in `src/`).
Holds the z-coordinates and their sort permutation, refreshed in place by
[`update_sorter!`](@ref) so repeated queries in an MD-style loop do not
reallocate. [`QuasiEwald.energy`](@ref)/[`force`](@ref)/[`force!`](@ref) on
[`QuasiEwaldLongPlan`](@ref) accept `z_list = sorter.z_list` (or compute
their own z-sort when no `z_list` is given at all).
"""
mutable struct ZSorter{T, TI}
    z_coords::Vector{T}
    z_list::Vector{TI}
end

function ZSorter(poses)
    z_coords = [p[3] for p in poses]
    return ZSorter(z_coords, sortperm(z_coords))
end

# The z-coordinates-only entry point the docstring above advertises. It used
# to be documented but not defined, so the call landed on `ZSorter(poses)`
# and raised `BoundsError: attempt to access Float64 at index [3]`. `copy` is
# not optional: `update_sorter!` writes into `sorter.z_coords` in place, so
# the sorter has to own its array rather than alias the caller's.
ZSorter(z_coords::Vector{T}) where {T<:Number} = ZSorter(copy(z_coords), sortperm(z_coords))

"Refresh `sorter` in place from the z-component of `poses` (no reallocation)."
function update_sorter!(sorter::ZSorter, poses)
    for i in eachindex(sorter.z_coords)
        sorter.z_coords[i] = poses[i][3]
    end
    sortperm!(sorter.z_list, sorter.z_coords)
    return sorter
end