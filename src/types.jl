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

struct QuasiEwaldShortInteraction{T, TI} <: ExTinyMD.AbstractInteraction
    # common used part
    γ_1::T
    γ_2::T
    ϵ_0::T
    L::NTuple{3, T}
    rbe::Bool
    accuracy::T
    α::T
    n_atoms::TI

    # short range part
    r_c::T
    n_t::TI
    gauss_para::GaussParameter{T}
end

QuasiEwaldShortInteraction(γ_1::T, γ_2::T, ϵ_0::T, L::NTuple{3, T}, rbe::Bool, accuracy::T, α::T, n_atoms::TI, r_c::T, n_t::TI) where {T<:Number, TI<:Integer} = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t, GaussParameter(n_t))

struct QuasiEwaldLongInteraction{T, TI} <: ExTinyMD.AbstractInteraction
    # common used part
    γ_1::T
    γ_2::T
    ϵ_0::T
    L::NTuple{3, T}
    rbe::Bool
    accuracy::T
    α::T
    n_atoms::TI

    # long range part
    k_c::T
    rbe_p::TI
    sum_k::T
    K_set::Vector{NTuple{3, T}}

    # divergent part
    k_0::T
    ringangles::RingAngles{T}

    # charge and coords
    q::Vector{T}
    mass::Vector{T}
    coords::Vector{Point{3, T}}
    acceleration::Vector{Point{3, T}}
end

function QuasiEwaldLongInteraction(γ_1::T, γ_2::T, ϵ_0::T, L::NTuple{3, T}, rbe::Bool, accuracy::T, α::T, n_atoms::TI, k_c::T, rbe_p::TI; Δk::T = π / sqrt(L[1] * L[2])) where{T<:Number, TI<:Integer}
    K_set, sum_k = rbe_sampling(L, α, accuracy)
    
    if γ_1 * γ_2 ≥ one(T)
        k_0 = log(γ_1 * γ_2) / (2 * L[3])
        ringangles = RingAngles(k_0, L[1], L[2], L[3], α, k_c, Δk)
    else
        k_0 = zero(T)
        ringangles = RingAngles(k_0)
    end

    q = zeros(T, n_atoms)
    mass = zeros(T, n_atoms)
    coords = Vector{Point{3, T}}(undef, n_atoms)
    acceleration = Vector{Point{3, T}}(undef, n_atoms)

    return QuasiEwaldLongInteraction{T, TI}(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, k_c, rbe_p, sum_k, K_set, k_0, ringangles, q, mass, coords, acceleration)
end

mutable struct SortingFinder{T, TI} <: ExTinyMD.AbstractNeighborFinder
    z_coords::Vector{T}
    z_list::Vector{TI}
end

function SortingFinder(info::SimulationInfo{T}) where {T<: Number}
    z_coords = [p_info.position[3] for p_info in info.particle_info]
    z_list = sortperm(z_coords)
    return SortingFinder{T, eltype(z_list)}(z_coords, z_list)
end

function ExTinyMD.update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: SortingFinder}
    n_atoms = length(neighborfinder.z_list)
    for i in 1:n_atoms
        neighborfinder.z_coords[i] = info.particle_info[i].position[3]
    end
    sortperm!(neighborfinder.z_list, neighborfinder.z_coords)
    return nothing
end