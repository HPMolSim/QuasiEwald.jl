export QuasiEwald_Fs!, QuasiEwald_Fs_pair, QuasiEwald_Fs_self

function Fsr_gauss_core(k::T, element::GreensElement{T}; l::Int = 0) where {T<:Number}
    f_sr_g = k * Gamma_1(k, element; l = l) * exp(- k*k / (4 * element.α)) * besselj1(k * element.ρ)
    return f_sr_g
end

function Fsr_point_core(k::T, element::GreensElement{T}; l::Int = 0) where {T<:Number}
    f_sr_p = k * Gamma_2(k, element; l = l) * besselj1(k * element.ρ)
    return f_sr_p
end

function Fsz_gauss_core(k::T, element::GreensElement{T}; l::Int = 0) where {T<:Number}
    f_sz_g = (exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) .* dz_Gamma_1(k, element; l = l)
    return f_sz_g
end

function Fsz_point_core(k::T, element::GreensElement{T}; l::Int = 0) where {T<:Number}
    f_sz_p = besselj0(k * element.ρ) .* dz_Gamma_2(k, element; l = l)
    return f_sz_p
end

function Fsz_self_gauss_core(k::T, element::GreensElement{T}; l::Int = 0) where {T<:Number}
    f_sz_g = (exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) .* dz_Gamma_self_1(k, element; l = l)
    return f_sz_g
end

function Fsz_self_point_core(k::T, element::GreensElement{T}; l::Int = 0) where {T<:Number}
    f_sz_p = besselj0(k * element.ρ) .* dz_Gamma_self_2(k, element; l = l)
    return f_sz_p
end

function QuasiEwald_Fs!(interaction::QuasiEwaldShortInteraction{T, TI}, neighborfinder::CellListDirQ2D{T, TI}, atoms::Vector{ExTinyMD.Atom{T}}, boundary::ExTinyMD.Boundary{T}, coords::Vector{Point{3, T}}, acceleration::Vector{Point{3, T}}) where {T<:Number, TI<:Integer}
    
    for (i, j, ρ) in neighborfinder.neighbor_list
        coord_1, coord_2, ρ_sq = position_checkQ2D(coords[i], coords[j], boundary, interaction.cutoff)
        if iszero(dist_sq)
            nothing
        else
            element = GreensElement(interaction.γ_1, interaction.γ_2, coord_1[3], coord_2[3], sqrt(ρ_sq), interaction.L[3], interaction.α, interaction.accuracy)
            q_1 = atoms[i].charge
            q_2 = atoms[j].charge
            force_i, force_j = QuasiEwald_Fs_pair(q_1, q_2, interaction.ϵ_0, element, coord_1, coord_2, interaction.gauss_para; single_mode = single_mode)
            acceleration[i] += force_i / atoms[i].mass
            acceleration[j] += force_j / atoms[j].mass
        end
    end

    for i in 1:interaction.n_atoms
        element = GreensElement(interaction.γ_1, interaction.γ_2, coords[i][3], interaction.L[3], interaction.α, interaction.accuracy)
        q = atoms[i].charge
        force_i = QuasiEwald_Fs_self(q, interaction.ϵ_0, element, interaction.gauss_para; single_mode = single_mode)

        acceleration[i] += force_i / atoms[i].mass
    end
    return nothing
end

function QuasiEwald_Fs_pair(q_1::T, q_2::T, ϵ_0::T, element::GreensElement{T}, coord_1::Point{3, T}, coord_2::Point{3, T}, gauss_para::GaussParameter{T}; single_mode::Bool = false) where T<:Number
    k_f1 = maximum(element.k_f1)
    k_f2 = maximum(element.k_f2)
    ρ = element.ρ

    # about the force in ρ direction
    Fsr_point_1 = Gauss_int(Fsr_point_core, gauss_para, element, region = (zero(T), k_f2))
    Fsr_point_2 = T(0.5) * sum(l -> element.b[l] * ρ / (element.a[l]^2 + ρ^2)^1.5, (1, 2, 3, 4))
    if single_mode == false
        Fsr_gauss = Gauss_int(Fsr_gauss_core, gauss_para, element, region = (zero(T), k_f1))
    else
        Fsr_gauss = zero(T)
    end

    Fsr = - Fsr_point_1 + Fsr_point_2 + Fsr_gauss
    Fsx = Fsr * (coord_1[1] - coord_2[1]) / ρ
    Fsy = Fsr * (coord_1[2] - coord_2[2]) / ρ
    
    # about the force in z direction
    Fsz_point_1 = Gauss_int_Tuple(Fsz_point_core, gauss_para, element, region = (zero(T), k_f2))
    
    a = element.a
    sa = element.sign_a
    b = element.b
    Fsz_point_2_temp = (
        b[1] * a[1] * sa[1] / (a[1]^2 + ρ^2)^1.5,
        b[2] * a[2] * sa[2] / (a[2]^2 + ρ^2)^1.5,
        b[3] * a[3] * sa[3] / (a[3]^2 + ρ^2)^1.5,
        b[4] * a[4] * sa[4] / (a[4]^2 + ρ^2)^1.5
    )
    Fsz_point_2 = (sum(Fsz_point_2_temp), dot((-one(T), one(T), one(T), -one(T)), Fsz_point_2_temp)) ./ T(2)
    if single_mode == false
        Fsz_gauss = Gauss_int_Tuple(Fsz_gauss_core, gauss_para, element, region = (zero(T), k_f1))
    else
        Fsz_gauss = (zero(T), zero(T))
    end
    Fsz = Fsz_point_1 .- Fsz_point_2 .- Fsz_gauss

    return (q_1 * q_2 / (2 * π * ϵ_0)) .* (Point(Fsx, Fsy, Fsz[1]), Point(-Fsx, -Fsy, Fsz[2]))
end

function QuasiEwald_Fs_self(q::T, ϵ_0::T, element::GreensElement{T}, gauss_para::GaussParameter{T}; single_mode::Bool = false) where T<:Number
    k_f1 = maximum(element.k_f1)
    k_f2 = maximum(element.k_f2)
    
    Fsz_point_1 = Gauss_int(Fsz_self_point_core, gauss_para, element, region = (zero(T), k_f2))
    a = element.a
    sa = element.sign_a
    b = element.b
    Fsz_point_2 = 0.5 * sum(l-> b[l] * sa[l] / a[l]^2, (2, 3))
    if single_mode == false
        Fsz_gauss = Gauss_int(Fsz_self_gauss_core, gauss_para, element, region = (zero(T), k_f1))
    else
        Fsz_gauss = zero(T)
    end

    Fsz = q^2 * Point(zero(T), zero(T), + Fsz_point_1 - Fsz_point_2 - Fsz_gauss) / (2 * π * ϵ_0)
    
    return Fsz
end