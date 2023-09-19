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
# now these two dz function will have tuple as return value
function dz_Gamma_1(k::T, element::GreensElement{T}; l::Int = 0) where T<:Number
    a = element.a
    sa = element.sign_a
    b = element.b

    dz_green_ui = (
            sa[1] * b[1] * exp(- k * a[1]), 
            sa[2] * b[2] * exp(- k * a[2]), 
            sa[3] * b[3] * exp(- k * a[3]), 
            sa[4] * b[4] * exp(- k * a[4])
            )
    dz_green_u1 = sum(dz_green_ui)
    dz_green_u2 = dot((-one(T), one(T), one(T), -one(T)), dz_green_ui)
    dz_green_u = (k/2) .* (dz_green_u1, dz_green_u2)

    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1

    dz_G_1 = dz_green_u ./ green_d
    return dz_G_1
end


function dz_Gamma_2(k::T, element::GreensElement{T}; l::Int = 0) where T<:Number
    return (element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z)) .* dz_Gamma_1(k, element)
end

function dz_Gamma_self_1(k::T, element::GreensElement{T}; l::Int = 0)::T where T<:Number
    a = element.a
    sa = element.sign_a
    b = element.b

    dz_green_u = (k/2) * (sa[2] * b[2] * exp(- k * a[2]) + sa[3] * b[3] * exp(- k * a[3]))

    green_d = element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z) - 1

    dz_G_1 = dz_green_u / green_d
    return dz_G_1
end


function dz_Gamma_self_2(k::T, element::GreensElement{T}; l::Int = 0)::T where T<:Number
    return (element.γ_1 * element.γ_2 * exp(- 2 * k * element.L_z)) * dz_Gamma_self_1(k, element)
end