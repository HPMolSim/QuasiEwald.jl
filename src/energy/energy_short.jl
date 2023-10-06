# the core functions are the integrands
function Es_gauss_core(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - T(1)

    if element.γ_1 * element.γ_2 ≤ 1
        E_s_g = Gamma_1(k, element) * exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ) / green_d
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        E_s_g = (Gamma_1(k, element) * exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ) - element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) * Gamma_1(k_0, element) * exp(- k_0*k_0 / (4 * element.α)) * besselj0(k_0 * element.ρ)) / green_d
    end

    return E_s_g
end

function Es_gauss_core(element::GreensElement{T}) where {T<:Number}
    if element.γ_1 * element.γ_2 ≤ 1
        E_s_g = zero(T)
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        E_s_g = Gamma_1(k_0, element) * exp(- k_0 * k_0 / (4 * element.α)) * besselj0(k_0 * element.ρ) * log(element.γ_1 * element.γ_2 - 1) / (2 * element.L_z)
    end
    return E_s_g
end

function Es_point_core(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - T(1)

    if element.γ_1 * element.γ_2 ≤ 1
        E_s_p = Gamma_2(k, element) * besselj0(k * element.ρ) / green_d
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        E_s_p = (Gamma_2(k, element) * besselj0(k * element.ρ) - element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) * Gamma_2(k_0, element) * besselj0(k_0 * element.ρ)) / green_d
    end
    
    return E_s_p
end

function Es_point_core(element::GreensElement{T}) where {T<:Number}
    if element.γ_1 * element.γ_2 ≤ 1
        E_s_p = zero(T)
    else
        k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
        E_s_p = Gamma_2(k_0, element) * besselj0(k_0 * element.ρ) * log(element.γ_1 * element.γ_2 - 1) / (2 * element.L_z)
    end

    return E_s_p
end

function QuasiEwald_Es(interaction::QuasiEwaldShortInteraction{T, TI}, neighbor::CellListDirQ2D{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}; single_mode::Bool = false) where {T<:Number, TI<:Integer}

    neighbor_list = neighbor.neighbor_list
    atoms = sys.atoms
    coords = info.coords

    energy_short = zero(T)

    for (i, j, ρ) in neighbor_list
        coord_1, coord_2, ρ_sq = position_checkQ2D(coords[i], coords[j], sys.boundary, interaction.r_c)
        if iszero(ρ_sq)
            nothing
        else
            element = GreensElement(interaction.γ_1, interaction.γ_2, coord_1[3], coord_2[3], sqrt(ρ_sq), interaction.L[3], interaction.α, interaction.accuracy)
            q_1 = atoms[i].charge
            q_2 = atoms[j].charge
            energy_short += QuaisEwald_Es_pair(q_1, q_2, interaction.ϵ_0, element, interaction.gauss_para; single_mode = single_mode)
        end
    end

    for i in 1:interaction.n_atoms
        element = GreensElement(interaction.γ_1, interaction.γ_2, coords[i][3], interaction.L[3], interaction.α, interaction.accuracy)
        q = atoms[i].charge
        energy_short += QuaisEwald_Es_self(q, interaction.ϵ_0, element, interaction.gauss_para; single_mode = single_mode)
    end

    return energy_short
end

function QuaisEwald_Es_pair(q_1::T, q_2::T, ϵ_0::T, element::GreensElement{T}, gauss_para::GaussParameter{T}; single_mode::Bool = false) where T<:Number
    k_f1 = maximum(element.k_f1)
    k_f2 = maximum(element.k_f2)
    Es_point_1 = Gauss_int(Es_point_core, gauss_para, element, region = (zero(T), k_f2)) + Es_point_core(element)
    Es_point_2 = 0.5 * sum(l -> element.b[l] / sqrt(element.a[l]^2 + element.ρ^2), (1, 2, 3, 4))

    if single_mode == false
        Es_gauss = Gauss_int(Es_gauss_core, gauss_para, element, region = (zero(T), k_f1)) + Es_gauss_core(element)
    else
        Es_gauss = zero(T)
    end

    Es_pair = q_1 * q_2 * (- Es_point_1 + Es_point_2 + Es_gauss) / (2π * ϵ_0)
    return Es_pair
end

function QuaisEwald_Es_self(q::T, ϵ_0::T, element::GreensElement{T}, gauss_para::GaussParameter{T}; single_mode::Bool = false) where T<:Number
    k_f1 = maximum(element.k_f1)
    k_f2 = maximum(element.k_f2)

    Es_point_1 = Gauss_int(Es_point_core, gauss_para, element, region = (zero(T), k_f2)) + Es_point_core(element)
    Es_point_2 = 0.5 * sum(l -> element.b[l] / element.a[l], (2, 3, 4))

    if single_mode == false
        Es_gauss = Gauss_int(Es_gauss_core, gauss_para, element, region = (zero(T), k_f1)) + Es_gauss_core(element)
    else
        Es_gauss = zero(T)
    end

    Es_self = q * q * (- Es_point_1 + Es_point_2 + Es_gauss) / (4π * ϵ_0)
    return Es_self
end