export Point, IcmSys, GaussParameter, GreensElement

struct Point{N,T}
    coo::NTuple{N,T}
end
Point(arg::T, args::T...) where T<:Number = Point((arg, args...))
Base.:(+)(x::Point{N,T}, y::Point{N,T}) where {N, T} = Point(x.coo .+ y.coo)
Base.:(-)(x::Point{N,T}, y::Point{N,T}) where {N, T} = Point(x.coo .- y.coo)
Base.:(-)(x::Point{N,T}) where {N, T} = Point(Base.:(-).(x.coo))
Base.adjoint(x::Point) = x
Base.:(*)(x::Number, y::Point) = Point(y.coo .* x)
Base.:(*)(y::Point, x::Number) = Point(y.coo .* x)
Base.iterate(x::Point, args...) = Base.iterate(x.coo, args...)
Base.getindex(x::Point, i::Int) = x.coo[i]

dist2(x::Number, y::Number) = abs2(x - y)
dist2(x::Point, y::Point) = sum(abs2, x - y)

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

    sign_a = (zero(T), -one(T), one(T), zero(T))
    k_f1 = sqrt.(4 * α^2 .* a.^2 .- 4 * α * log(accuracy)) .- 2 * α .* a
    k_f2 = - log(accuracy) ./ (2 * L_z .+ a)

    return GreensElement{T}(γ_1, γ_2, ρ, a, b, sign_a, L_z, α, k_f1, k_f2)
end