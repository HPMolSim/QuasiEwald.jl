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
# ============================================================================
# Framework-free core query (Task 3).
# ============================================================================

"""
    QuasiEwald.energy(plan::QuasiEwaldShortPlan, poses, charges; neighbor_list = nothing) -> T

Short-range (real-space) energy from plain array-of-structs positions and
charges -- no ExTinyMD type constructed, neither argument mutated.

Pass `neighbor_list` (an iterable of `(i, j, ...)` candidate pairs, e.g. a
`CellListMap` neighbor list or the ExTinyMD extension's `CellListQ2D`) to
reuse one already maintained elsewhere; each pair is treated as a candidate
only -- the in-plane distance is always recomputed here via the true
minimum image (`_min_image_q2d`), so a supplied list's own reported
distance is ignored, matching the discipline ExTinyMD's stdlib uses for the
same reason (a supplied finder's metric may not match this plan's own,
e.g. an in-plane-only cell list reporting `r = 0`).

With no `neighbor_list`, every pair is tested directly (`O(n_atoms^2)`) --
this plan owns no persistent cell list of its own; see the package README
for why that trade was made.
"""
function energy(plan::QuasiEwaldShortPlan{T}, poses, charges; neighbor_list = nothing) where {T}
    n_atoms = plan.n_atoms
    energy_short = zero(T)
    r_c_sq = plan.r_c^2

    if neighbor_list === nothing
        for i in 1:n_atoms, j in (i + 1):n_atoms
            energy_short += _short_pair_energy(plan, poses, charges, i, j, r_c_sq)
        end
    else
        for pair in neighbor_list
            i, j = pair[1], pair[2]
            energy_short += _short_pair_energy(plan, poses, charges, i, j, r_c_sq)
        end
    end

    for i in 1:n_atoms
        element = GreensElement(plan.γ_1, plan.γ_2, poses[i][3], plan.L[3], plan.α, plan.accuracy)
        energy_short += QuaisEwald_Es_self(charges[i], plan.ϵ_0, element, plan.gauss_para)
    end

    return energy_short
end

"Candidate-pair energy contribution, zero unless within `r_c_sq` after the true minimum-image correction."
@inline function _short_pair_energy(plan::QuasiEwaldShortPlan{T}, poses, charges, i, j, r_c_sq::T) where {T}
    coord_1, coord_2, ρ_sq = _min_image_q2d(poses[i], poses[j], plan.L)
    if ρ_sq ≥ r_c_sq
        return zero(T)
    end
    element = GreensElement(plan.γ_1, plan.γ_2, coord_1[3], coord_2[3], sqrt(ρ_sq), plan.L[3], plan.α, plan.accuracy)
    return QuaisEwald_Es_pair(charges[i], charges[j], plan.ϵ_0, element, plan.gauss_para)
end
