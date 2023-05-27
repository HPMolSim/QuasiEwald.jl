export Gamma_1, Gamma_2, dz_Gamma_1, dz_Gamma_2

function Gamma_1(k::T, element::GreensElement{T}; l::Int = 0) where T<:Number
    if l == 0
        green_u = sum(i -> element.b[i] * exp(- k * element.a[i]), (1, 2, 3, 4)) / 2
    else
        green_u = element.b[l] * exp(- k * element.a[l]) / 2
    end
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    G_1 = green_u / green_d
    return G_1
end

function Gamma_2(k::T, element::GreensElement{T}; l::Int = 0) where T<:Number
    if l == 0
        green_u = sum(i -> element.b[i] * exp(- k * element.a[i]), (1, 2 ,3, 4)) / 2
    else
        green_u = element.b[l] * exp(- k * element.a[l]) / 2
    end
    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    G_2 = element.γ_1 * element.γ_2 * green_u * exp(- 2.0 * k * element.L_z) / green_d
    return G_2
end


# these functions define the dz_Gamma1/2, which are used to calculate the forces
function dz_Gamma_1(k::T, element::GreensElement{T}; l::Int = 0) where T<:Number
    if l == 0
        dz_green_u = sum(i -> k * element.sign_a[i] * element.b[i] * exp(- k * element.a[i]), (1, 2, 3, 4)) / 2
    else
        dz_green_u = k * element.sign_a[l] * element.b[l] * exp(- k * element.a[l]) / 2
    end

    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1

    dz_G_1 = dz_green_u / green_d
    return dz_G_1
end

function dz_Gamma_2(k::T, element::GreensElement{T}; l::Int = 0) where T<:Number
    if l == 0
        dz_green_u = sum(i -> k * element.sign_a[i] * element.b[i] * exp(- k * element.a[i]), (1, 2, 3, 4)) / 2
    else
        dz_green_u = k * element.sign_a[l] * element.b[l] * exp(- k * element.a[l]) / 2
    end

    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1
    dz_G_2 = element.γ_1 * element.γ_2 * dz_green_u * exp(- 2 * k * element.L_z) / green_d

    return dz_G_2
end