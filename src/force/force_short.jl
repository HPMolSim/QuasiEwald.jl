function Fsr_gauss_core(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - T(1)

    if element.γ_1 * element.γ_2 ≤ 1
        f_sr_g = k * Gamma_1(k, element) * exp(- k*k / (4 * element.α)) * besselj1(k * element.ρ) / green_d
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sr_g = (k * Gamma_1(k, element) * exp(- k*k / (4 * element.α)) * besselj1(k * element.ρ) - element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) * k_0 * Gamma_1(k_0, element) * exp(- k_0*k_0 / (4 * element.α)) * besselj1(k_0 * element.ρ)) / green_d
    end
    return f_sr_g
end

function Fsr_gauss_core(element::GreensElement{T}) where {T<:Number}
    if element.γ_1 * element.γ_2 ≤ 1
        f_sr_g = zero(T)
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sr_g = k_0 * Gamma_1(k_0, element) * exp(- k_0*k_0 / (4 * element.α)) * besselj1(k_0 * element.ρ) * log(element.γ_1 * element.γ_2 - 1) / (2 * element.L_z)
    end
    return f_sr_g
end

function Fsr_point_core(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - T(1)

    if element.γ_1 * element.γ_2 ≤ 1
        f_sr_p = k * Gamma_2(k, element) * besselj1(k * element.ρ) / green_d
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sr_p = (k * Gamma_2(k, element) * besselj1(k * element.ρ) - element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) * k_0 * Gamma_2(k_0, element) * besselj1(k_0 * element.ρ)) / green_d
    end

    return f_sr_p
end

function Fsr_point_core(element::GreensElement{T}) where {T<:Number}

    if element.γ_1 * element.γ_2 ≤ 1
        f_sr_p = zero(T)
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sr_p = k_0 * Gamma_2(k_0, element) * besselj1(k_0 * element.ρ) * log(element.γ_1 * element.γ_2 - 1) / (2 * element.L_z)
    end
    return f_sr_p
end

function Fsz_gauss_core(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1

    if element.γ_1 * element.γ_2 ≤ 1
        f_sz_g = (exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) .* dz_Gamma_1(k, element) ./ green_d
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sz_g = ((exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) .* dz_Gamma_1(k, element) .- (element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) * (exp(- k_0*k_0 / (4 * element.α)) * besselj0(k_0 * element.ρ))) .* dz_Gamma_1(k_0, element)) ./ green_d
    end

    return f_sz_g
end

function Fsz_gauss_core(element::GreensElement{T}) where {T<:Number}
    if element.γ_1 * element.γ_2 ≤ 1
        f_sz_g = zero(T)
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sz_g = ((exp(- k_0*k_0 / (4 * element.α)) * besselj0(k_0 * element.ρ)) * log(element.γ_1 * element.γ_2 - 1) / (2 * element.L_z)) .* dz_Gamma_1(k_0, element)
    end
    return f_sz_g
end

function Fsz_point_core(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    if element.γ_1 * element.γ_2 ≤ 1
        f_sz_p = besselj0(k * element.ρ) .* dz_Gamma_2(k, element) ./ green_d
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sz_p = (besselj0(k * element.ρ) .* dz_Gamma_2(k, element) .- (element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z)) * besselj0(k_0 * element.ρ) .* dz_Gamma_2(k_0, element) ) ./ green_d
    end
    return f_sz_p
end

function Fsz_point_core(element::GreensElement{T}) where {T<:Number}
    if element.γ_1 * element.γ_2 ≤ 1
        f_sz_p = zero(T)
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sz_p = (besselj0(k_0 * element.ρ) * log(element.γ_1 * element.γ_2 - 1) / (2 * element.L_z)) .* dz_Gamma_2(k_0, element)
    end
    return f_sz_p
end

function Fsz_self_gauss_core(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    if element.γ_1 * element.γ_2 ≤ 1
        f_sz_g = (exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) * dz_Gamma_self_1(k, element) / green_d
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sz_g = ((exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) * dz_Gamma_self_1(k, element) - element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) * (exp(- k_0*k_0 / (4 * element.α)) * besselj0(k_0 * element.ρ)) * dz_Gamma_self_1(k_0, element)) / green_d
    end
    return f_sz_g
end

function Fsz_self_gauss_core(element::GreensElement{T}) where {T<:Number}
    if element.γ_1 * element.γ_2 ≤ 1
        f_sz_g = zero(T)
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sz_g = (exp(- k_0*k_0 / (4 * element.α)) * besselj0(k_0 * element.ρ)) * dz_Gamma_self_1(k_0, element) * log(element.γ_1 * element.γ_2 - 1) / (2 * element.L_z)
    end
    return f_sz_g
end

function Fsz_self_point_core(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    if element.γ_1 * element.γ_2 ≤ 1
        f_sz_p = besselj0(k * element.ρ) .* dz_Gamma_self_2(k, element) / green_d
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sz_p = (besselj0(k * element.ρ) .* dz_Gamma_self_2(k, element) - element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) * besselj0(k_0 * element.ρ) .* dz_Gamma_self_2(k_0, element)) / green_d
    end
    return f_sz_p
end

function Fsz_self_point_core( element::GreensElement{T}) where {T<:Number}
    if element.γ_1 * element.γ_2 ≤ 1
        f_sz_p = zero(T)
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        f_sz_p = besselj0(k_0 * element.ρ) * dz_Gamma_self_2(k_0, element) * log(element.γ_1 * element.γ_2 - 1) / (2 * element.L_z)
    end
    return f_sz_p
end

function QuasiEwald_Fs!(interaction::QuasiEwaldShortInteraction{T, TI}, neighborfinder::CellListQ2D{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}

    atoms = sys.atoms
    
    for (i, j, ρ) in neighborfinder.neighbor_list
        id_i = info.particle_info[i].id
        id_j = info.particle_info[j].id
        coord_1, coord_2, ρ_sq = position_checkQ2D(info.particle_info[i].position, info.particle_info[j].position, sys.boundary, interaction.r_c)
        if iszero(ρ_sq)
            nothing
        else
            element = GreensElement(interaction.γ_1, interaction.γ_2, coord_1[3], coord_2[3], sqrt(ρ_sq), interaction.L[3], interaction.α, interaction.accuracy)
            q_1 = atoms[id_i].charge
            q_2 = atoms[id_j].charge
            force_i, force_j = QuasiEwald_Fs_pair(q_1, q_2, interaction.ϵ_0, element, coord_1, coord_2, interaction.gauss_para)
            info.particle_info[i].acceleration += force_i / atoms[id_i].mass
            info.particle_info[j].acceleration += force_j / atoms[id_j].mass
        end
    end

    for p_info in info.particle_info
        id_i = p_info.id
        element = GreensElement(interaction.γ_1, interaction.γ_2, p_info.position[3], interaction.L[3], interaction.α, interaction.accuracy)
        q = atoms[p_info.id].charge
        force_i = QuasiEwald_Fs_self(q, interaction.ϵ_0, element, interaction.gauss_para)

        p_info.acceleration += force_i / atoms[id_i].mass
    end
    
    return nothing
end

function QuasiEwald_Fs_pair(q_1::T, q_2::T, ϵ_0::T, element::GreensElement{T}, coord_1::Point{3, T}, coord_2::Point{3, T}, gauss_para::GaussParameter{T}; single_mode::Bool = false) where T<:Number
    k_f1 = maximum(element.k_f1)
    k_f2 = maximum(element.k_f2)
    ρ = element.ρ

    # about the force in ρ direction
    Fsr_point_1 = Gauss_int(Fsr_point_core, gauss_para, element, region = (zero(T), k_f2)) + Fsr_point_core(element)
    Fsr_point_2 = T(0.5) * sum(l -> element.b[l] * ρ / (element.a[l]^2 + ρ^2)^1.5, (1, 2, 3, 4))
    if single_mode == false
        Fsr_gauss = Gauss_int(Fsr_gauss_core, gauss_para, element, region = (zero(T), k_f1)) + Fsr_gauss_core(element)
    else
        Fsr_gauss = zero(T)
    end

    Fsr = - Fsr_point_1 + Fsr_point_2 + Fsr_gauss
    Fsx = Fsr * (coord_1[1] - coord_2[1]) / ρ
    Fsy = Fsr * (coord_1[2] - coord_2[2]) / ρ
    
    # about the force in z direction
    Fsz_point_1 = Gauss_int_Tuple(Fsz_point_core, gauss_para, element, region = (zero(T), k_f2)) .+ Fsz_point_core(element)
    
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
        Fsz_gauss = Gauss_int_Tuple(Fsz_gauss_core, gauss_para, element, region = (zero(T), k_f1)) .+ Fsz_gauss_core(element)
    else
        Fsz_gauss = (zero(T), zero(T))
    end
    Fsz = Fsz_point_1 .- Fsz_point_2 .- Fsz_gauss

    return (q_1 * q_2 / (2 * π * ϵ_0)) .* (Point(Fsx, Fsy, Fsz[1]), Point(-Fsx, -Fsy, Fsz[2]))

end

function QuasiEwald_Fs_self(q::T, ϵ_0::T, element::GreensElement{T}, gauss_para::GaussParameter{T}; single_mode::Bool = false) where T<:Number
    k_f1 = maximum(element.k_f1)
    k_f2 = maximum(element.k_f2)
    
    Fsz_point_1 = Gauss_int(Fsz_self_point_core, gauss_para, element, region = (zero(T), k_f2)) + Fsz_self_point_core(element)
    a = element.a
    sa = element.sign_a
    b = element.b
    Fsz_point_2 = 0.5 * sum(l-> b[l] * sa[l] / a[l]^2, (2, 3))
    if single_mode == false
        Fsz_gauss = Gauss_int(Fsz_self_gauss_core, gauss_para, element, region = (zero(T), k_f1)) + Fsz_self_gauss_core(element)
    else
        Fsz_gauss = zero(T)
    end

    Fsz = q^2 * Point(zero(T), zero(T), + Fsz_point_1 - Fsz_point_2 - Fsz_gauss) / (2 * π * ϵ_0)
    
    return Fsz
end