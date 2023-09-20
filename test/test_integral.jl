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


@testset "Gauss Integrator for divergent integrands (energy)" begin
    γ_1 = 2.0
    γ_2 = 3.0
    L_z = 10.0
    z_i = 1.1
    z_j = 2.2
    ρ = 3.3
    accuracy = 1e-8
    α = 0.9

    element = GreensElement(γ_1, γ_2, z_i, z_j, ρ, L_z, α, accuracy)
    gauss_para_pv = GaussParameter(50)
    gauss_para_direct = GaussParameter(200)

    k_f1 = maximum(element.k_f1)
    k_f2 = maximum(element.k_f2)
    k_0 = log(element.γ_1 * element.γ_2) / (2 * element.L_z)

    energy_gauss_pv = Gauss_int(Es_gauss_core, gauss_para_pv, element, region = (0.0, k_f1)) + Es_gauss_core(element)

    energy_gauss_direct = Gauss_int(Es_gauss_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) + Gauss_int(Es_gauss_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f1))

    energy_point_pv = Gauss_int(Es_point_core, gauss_para_pv, element, region = (0.0, k_f2)) + Es_point_core(element)

    energy_point_direct = Gauss_int(Es_point_core_direct, gauss_para_direct, element, region = (0.0, 2 * k_0)) + Gauss_int(Es_point_core_direct, gauss_para_direct, element, region = (2 * k_0, k_f2))

    @test isapprox(energy_gauss_pv, energy_gauss_direct, atol = 1e-6)
    @test isapprox(energy_point_pv, energy_point_direct, atol = 1e-6)
end