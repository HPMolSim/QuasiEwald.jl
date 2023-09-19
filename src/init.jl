export RBE_α, QuasiEwaldRbeInit

function RBE_α(n_atoms::TI, L_x::T, L_y::T, n_t::TI, accuracy::T, rbe_p::TI) where {T<:Number, TI<:Integer}
    return ((n_atoms / L_x / L_y) * 26 * n_t * π / (12 * rbe_p * (accuracy)^(2/3)))^(1.5)
end

function QuasiEwaldRbeInit(info::SimulationInfo{T}, n_atoms::TI, L_x::T, L_y::T, n_t::TI, accuracy::T, rbe_p::TI, γ_1::T, γ_2::T, ϵ_0::T) where{T<:Number, TI<:Integer}
    coords = info.coords

    α = RBE_α(n_atoms, L_x, L_y, n_t, accuracy, rbe_p)
    k_c = sqrt(- 4 * α * log(accuracy))
    r_c = (α * accuracy)^(-1/3)

    intershort = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, r_c, n_t)
    findershort = CellListDirQ2D(info, r_c + 1.0, boundary, 100)
    interlong = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, k_c, rbe_p)
    finderlong = SortingFinder(coords)

    return intershort, findershort, interlong, finderlong
end