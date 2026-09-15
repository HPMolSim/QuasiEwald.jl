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


function QuasiEwald_Fs_pair(q_1::T, q_2::T, ϵ_0::T, element::GreensElement{T}, coord_1, coord_2, gauss_para::GaussParameter{T}; single_mode::Bool = false) where T<:Number
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
    # `Fsr` is the RADIAL magnitude; the in-plane components need the unit
    # vector (dx, dy)/ρ. At ρ == 0 that is 0/0: `Fsr` itself vanishes there
    # (every term carries either `besselj1(k*ρ)` or an explicit factor of ρ),
    # so the quotient is NaN, not a genuine singularity. The limit is zero and
    # approached linearly -- measured Fx = 1.088e-5, 1.088e-7, 1.088e-10,
    # 1.088e-13 at dx = 1e-1, 1e-3, 1e-6, 1e-9 -- which is also what symmetry
    # demands, since dx and dy are identically zero. Fsz is finite and
    # continuous through ρ = 0 and is left alone.
    #
    # ρ == 0 is reachable and not exotic: any two particles sharing an (x, y)
    # column, and any pair whose in-plane separation is an exact multiple of
    # Lx or Ly (so the minimum image wraps to zero) -- i.e. any lattice or
    # grid initialisation. Before this package was decoupled, ExTinyMD's
    # `position_checkQ2D` returned its all-zero sentinel for such a pair and
    # every caller's `iszero(ρ_sq)` guard skipped it, which masked this by
    # dropping the pair's real z-force too. The explicit `ρ_sq ≥ r_c^2` test
    # that replaced the sentinel keeps the pair, so the guard has to be here.
    if iszero(ρ)
        Fsx = zero(T)
        Fsy = zero(T)
    else
        Fsx = Fsr * (coord_1[1] - coord_2[1]) / ρ
        Fsy = Fsr * (coord_1[2] - coord_2[2]) / ρ
    end
    
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

    return (q_1 * q_2 / (2 * π * ϵ_0)) .* (SVector{3, T}(Fsx, Fsy, Fsz[1]), SVector{3, T}(-Fsx, -Fsy, Fsz[2]))

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

    Fsz = q^2 * SVector{3, T}(zero(T), zero(T), + Fsz_point_1 - Fsz_point_2 - Fsz_gauss) / (2 * π * ϵ_0)
    
    return Fsz
end
# ============================================================================
# Framework-free core queries (Task 3).
# ============================================================================

"Candidate-pair force contribution on `i` and `j`; both zero unless within `r_c_sq` after the true minimum-image correction."
@inline function _short_pair_force(plan::QuasiEwaldShortPlan{T}, poses, charges, i, j, r_c_sq::T) where {T}
    coord_1, coord_2, ρ_sq = _min_image_q2d(poses[i], poses[j], plan.L)
    if ρ_sq ≥ r_c_sq
        z = SVector{3, T}(zero(T), zero(T), zero(T))
        return z, z
    end
    element = GreensElement(plan.γ_1, plan.γ_2, coord_1[3], coord_2[3], sqrt(ρ_sq), plan.L[3], plan.α, plan.accuracy)
    return QuasiEwald_Fs_pair(charges[i], charges[j], plan.ϵ_0, element, coord_1, coord_2, plan.gauss_para)
end

"""
    QuasiEwald.force!(F, plan::QuasiEwaldShortPlan, poses, charges; neighbor_list = nothing) -> F
    QuasiEwald.force(plan::QuasiEwaldShortPlan, poses, charges; neighbor_list = nothing) -> Vector{SVector{3,T}}

Short-range (real-space) force from plain array-of-structs positions and
charges, written into `F` (filled, not accumulated into). Neither `poses`
nor `charges` is mutated. See [`QuasiEwald.energy`](@ref) for the
`neighbor_list` contract (candidate pairs only, distance always recomputed
here) and for why no `neighbor_list` means an O(n_atoms^2) direct pair loop.
"""
function force!(F, plan::QuasiEwaldShortPlan{T}, poses, charges; neighbor_list = nothing) where {T}
    n_atoms = plan.n_atoms
    r_c_sq = plan.r_c^2
    fill!(F, SVector{3, T}(zero(T), zero(T), zero(T)))

    if neighbor_list === nothing
        for i in 1:n_atoms, j in (i + 1):n_atoms
            force_i, force_j = _short_pair_force(plan, poses, charges, i, j, r_c_sq)
            F[i] += force_i
            F[j] += force_j
        end
    else
        for pair in neighbor_list
            i, j = pair[1], pair[2]
            force_i, force_j = _short_pair_force(plan, poses, charges, i, j, r_c_sq)
            F[i] += force_i
            F[j] += force_j
        end
    end

    for i in 1:n_atoms
        element = GreensElement(plan.γ_1, plan.γ_2, poses[i][3], plan.L[3], plan.α, plan.accuracy)
        F[i] += QuasiEwald_Fs_self(charges[i], plan.ϵ_0, element, plan.gauss_para)
    end

    return F
end

"Allocating form of [`QuasiEwald.force!`](@ref)."
function force(plan::QuasiEwaldShortPlan{T}, poses, charges; kwargs...) where {T}
    F = [SVector{3, T}(zero(T), zero(T), zero(T)) for _ in 1:plan.n_atoms]
    return force!(F, plan, poses, charges; kwargs...)
end
