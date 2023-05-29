export IcmSys, GaussParameter, GreensElement, QuasiEwaldShortInteraction, QuasiEwaldLongInteraction

struct IcmSys{T, R}
    γ::NTuple{2, T} # (γ_up, γ_down)
    L::NTuple{3, T} # (Lx, Ly, Lz)
    N_real::R
    N_img::R
end

IcmSys(γ::NTuple{2, Float64}, L::NTuple{3, Float64}, N_real::Int, N_img::Int) = IcmSys{Float64, Int}(γ, L, N_real, N_img)

struct GaussParameter{N, T}
    sw::SVector{N, NTuple{2, T}}
end

GaussParameter(Step::Int) = GaussParameter{Float64}(SVector{Step, NTuple{2, Float64}}([tuple(legendre(Step)[1][i], legendre(Step)[2][i]) for i in 1:Step]))

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

    sign_a = (zero(T), -one(T), one(T), zero(T))
    k_f1 = sqrt.(4 * α^2 .* a.^2 .- 4 * α * log(accuracy)) .- 2 * α .* a
    k_f2 = - log(accuracy) ./ (2 * L_z .+ a)

    return GreensElement{T}(γ_1, γ_2, ρ, a, b, sign_a, L_z, α, k_f1, k_f2)
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
    # Prob::ProbabilityWeights{T, T, Vector{T}}
end

function QuasiEwaldLongInteraction(γ_1::T, γ_2::T, ϵ_0::T, L::NTuple{3, T}, rbe::Bool, accuracy::T, α::T, n_atoms::TI, k_c::T, rbe_p::TI) where{T<:Number, TI<:Integer}
    K_set, sum_k = rbe_sampling(L, α, accuracy)

    return QuasiEwaldLongInteraction{T, TI}(γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, k_c, rbe_p, sum_k, K_set, Prob)
end

mutable struct SortingFinder{T, TI} <: ExTinyMD.AbstractNeighborFinder
    z_coords::Vector{T}
    z_list::Vector{TI}
end

function SortingFinder(coords::Vector{Point{3, T}}) where {T<: Number}
    z_coords = [coord[3] for coord in coords]
    z_list = sortperm(z_coords)
    return SortingFinder{T, eltype(z_list)}(z_coords, z_list)
end

function ExTinyMD.update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: SortingFinder}
    n_atoms = length(neighborfinder.z_list)
    for i in 1:n_atoms
        neighborfinder.z_coords[i] = info.coords[i][3]
    end
    sortperm!(neighborfinder.z_list, neighborfinder.z_coords)
    return nothing
end