export Gauss_int

@inline function Gauss_int(integrand::Function, Gaussian::GaussParameter, para::GreensElement; region::NTuple{2, T} = (-1.0, 1.0)) where {T <: Number}
    
    a, b = region
    result = sum(i->integrand((b + a) / 2 + (b - a) * i[1] / 2, para; l = 0) * (b - a) * i[2] / 2, Gaussian.sw; init=0.0)

    return result
end