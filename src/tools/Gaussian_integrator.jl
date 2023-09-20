export Gauss_int, Gauss_int_Tuple

@inline function Gauss_int(integrand::Function, Gaussian::GaussParameter{T}, para::GreensElement{T}; region::NTuple{2, T} = (-1.0, 1.0)) where {T <: Number}

    a, b = region
    result = sum(i->integrand((b + a) / 2 + (b - a) * i[1] / 2, para) * (b - a) * i[2] / 2, Gaussian.sw; init= zero(T))
    
    return result
end

@inline function Gauss_int_Tuple(integrand::Function, Gaussian::GaussParameter{T}, para::GreensElement{T}; region::NTuple{2, T} = (-1.0, 1.0)) where {T <: Number}
    
    a, b = region
    result_1 = zero(T)
    result_2 = zero(T)
    temp_1 = zero(T)
    temp_2 = zero(T)
    for i in Gaussian.sw
        temp_1, temp_2 = integrand((b + a) / T(2) + (b - a) * i[1] / T(2), para)
        temp = (b - a) * i[2] / T(2)
        result_1 += temp_1 * temp
        result_2 += temp_2 * temp
    end

    return (result_1, result_2)
end