function Es_gauss_core_direct(k::T, element::GreensElement{T}) where {T<:Number}

    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - T(1)
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    if k ≤ 2 * k_0
        E_s_g = Gamma_1(k, element) * exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ) / green_d + Gamma_1(k_0, element) * exp(- k_0*k_0 / (4 * element.α)) * besselj0(k_0 * element.ρ) / (2 * element.L_z * (k - k_0))
    else
        E_s_g = Gamma_1(k, element) * exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ) / green_d
    end

    return E_s_g
end

function Es_point_core_direct(k::T, element::GreensElement{T}) where {T<:Number}

    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - T(1)
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    if k ≤ 2 * k_0
        E_s_p = Gamma_2(k, element) * besselj0(k * element.ρ) / green_d + Gamma_2(k_0, element) * besselj0(k_0 * element.ρ) / (2 * element.L_z * (k - k_0))
    else
        E_s_p = Gamma_2(k, element) * besselj0(k * element.ρ) / green_d
    end

    return E_s_p
end

function Fsr_gauss_core_direct(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - T(1)
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    if k ≤ 2 * k_0
        f_sr_g = k * Gamma_1(k, element) * exp(- k*k / (4 * element.α)) * besselj1(k * element.ρ) / green_d + k_0 * Gamma_1(k_0, element) * exp(- k_0*k_0 / (4 * element.α)) * besselj1(k_0 * element.ρ) / (2 * element.L_z * (k - k_0))
    else
        f_sr_g = k * Gamma_1(k, element) * exp(- k*k / (4 * element.α)) * besselj1(k * element.ρ) / green_d
    end

    return f_sr_g
end

function Fsr_point_core_direct(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - T(1)
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    if k ≤ 2 * k_0
        f_sr_p = k * Gamma_2(k, element) * besselj1(k * element.ρ) / green_d + k_0 * Gamma_2(k_0, element) * besselj1(k_0 * element.ρ) / (2 * element.L_z * (k - k_0))
    else
        f_sr_p = k * Gamma_2(k, element) * besselj1(k * element.ρ) / green_d
    end

    return f_sr_p
end

function Fsz_gauss_core_direct(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    if k ≤ 2 * k_0
        f_sz_g = (exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) .* dz_Gamma_1(k, element) ./ green_d .+ (exp(- k_0^2 / (4 * element.α)) * besselj0(k_0 * element.ρ)) .* dz_Gamma_1(k_0, element) ./ (2 * element.L_z * (k - k_0))
    else
        f_sz_g = (exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) .* dz_Gamma_1(k, element) ./ green_d
    end

    return f_sz_g
end

function Fsz_point_core_direct(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)
    
    if k ≤ 2 * k_0
        f_sz_p = besselj0(k * element.ρ) .* dz_Gamma_2(k, element) ./ green_d .+ besselj0(k_0 * element.ρ) .* dz_Gamma_2(k_0, element) ./ (2 * element.L_z * (k - k_0))
    else
        f_sz_p = besselj0(k * element.ρ) .* dz_Gamma_2(k, element) ./ green_d
    end

    return f_sz_p
end

function Fsz_self_gauss_core_direct(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    f_szs_g = (exp(- k*k / (4 * element.α)) * besselj0(k * element.ρ)) * dz_Gamma_self_1(k, element) / green_d

    if k ≤ 2 * k_0
        f_szs_g += (exp(- k_0*k_0 / (4 * element.α)) * besselj0(k_0 * element.ρ)) * dz_Gamma_self_1(k_0, element) / (2 * element.L_z * (k - k_0))
    end

    return f_szs_g
end

function Fsz_self_point_core_direct(k::T, element::GreensElement{T}) where {T<:Number}
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    f_szs_p = besselj0(k * element.ρ) .* dz_Gamma_self_2(k, element) / green_d

    if k ≤ 2 * k_0
        f_szs_p += besselj0(k_0 * element.ρ) .* dz_Gamma_self_2(k_0, element) / (2 * element.L_z * (k - k_0))
    end

    return f_szs_p
end


@testset "P.V. Integrator" begin
    γ_1 = 2.0
    γ_2 = 3.0
    L_z = 10.0
    z_i = 1.1
    z_j = 2.2
    ρ = 3.3
    accuracy = 1e-8
    α = 0.9

    element = GreensElement(γ_1, γ_2, z_i, z_j, ρ, L_z, α, accuracy)
    gauss_para_pv = GaussParameter(30)
    gauss_para_direct = GaussParameter(100)

    k_f1 = maximum(element.k_f1)
    k_f2 = maximum(element.k_f2)
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    energy_gauss_pv = Gauss_int(Es_gauss_core, gauss_para_pv, element, region = (0.0, k_f1)) + Es_gauss_core(element)
    energy_gauss_direct = Gauss_int(Es_gauss_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) + Gauss_int(Es_gauss_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f1))

    energy_point_pv = Gauss_int(Es_point_core, gauss_para_pv, element, region = (0.0, k_f2)) + Es_point_core(element)
    energy_point_direct = Gauss_int(Es_point_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) + Gauss_int(Es_point_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f2))

    fsr_gauss_pv = Gauss_int(Fsr_gauss_core, gauss_para_pv, element, region = (0.0, k_f1)) + Fsr_gauss_core(element)
    fsr_gauss_direct = Gauss_int(Fsr_gauss_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) + Gauss_int(Fsr_gauss_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f1))

    fsr_point_pv = Gauss_int(Fsr_point_core, gauss_para_pv, element, region = (0.0, k_f2)) + Fsr_point_core(element)
    fsr_point_direct = Gauss_int(Fsr_point_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) + Gauss_int(Fsr_point_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f2))

    fsz_gauss_pv = Gauss_int_Tuple(Fsz_gauss_core, gauss_para_pv, element, region = (0.0, k_f1)) .+ Fsz_gauss_core(element)
    fsz_gauss_direct = Gauss_int_Tuple(Fsz_gauss_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) .+ Gauss_int_Tuple(Fsz_gauss_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f1))

    fsz_point_pv = Gauss_int_Tuple(Fsz_point_core, gauss_para_pv, element, region = (0.0, k_f2)) .+ Fsz_point_core(element)
    fsz_point_direct = Gauss_int_Tuple(Fsz_point_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) .+ Gauss_int_Tuple(Fsz_point_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f2))

    fszs_gauss_pv = Gauss_int(Fsz_self_gauss_core, gauss_para_pv, element, region = (0.0, k_f1)) + Fsz_self_gauss_core(element)
    fszs_gauss_direct = Gauss_int(Fsz_self_gauss_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) + Gauss_int(Fsz_self_gauss_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f1))

    fszs_point_pv = Gauss_int(Fsz_self_point_core, gauss_para_pv, element, region = (0.0, k_f2)) + Fsz_self_point_core(element)
    fszs_point_direct = Gauss_int(Fsz_self_point_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) + Gauss_int(Fsz_self_point_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f2))


    @test isapprox(energy_gauss_pv, energy_gauss_direct, atol = 1e-6)
    @test isapprox(energy_point_pv, energy_point_direct, atol = 1e-6)
    @test isapprox(fsr_gauss_pv, fsr_gauss_direct, atol = 1e-6)
    @test isapprox(fsr_point_pv, fsr_point_direct, atol = 1e-6)
    @test isapprox(fsz_gauss_pv[1], fsz_gauss_direct[1], atol = 1e-6)
    @test isapprox(fsz_gauss_pv[2], fsz_gauss_direct[2], atol = 1e-6)
    @test isapprox(fsz_point_pv[1], fsz_point_direct[1], atol = 1e-6)
    @test isapprox(fsz_point_pv[2], fsz_point_direct[2], atol = 1e-6)
    @test isapprox(fszs_gauss_pv, fszs_gauss_direct, atol = 1e-6)
    @test isapprox(fszs_point_pv, fszs_point_direct, atol = 1e-6)
end