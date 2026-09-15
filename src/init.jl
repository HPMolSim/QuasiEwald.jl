function RBE_α(n_atoms::TI, L_x::T, L_y::T, n_t::TI, accuracy::T, rbe_p::TI) where {T<:Number, TI<:Integer}
    return ((n_atoms / L_x / L_y) * 26 * n_t * π / (12 * rbe_p * (accuracy)^(2/3)))^(1.5)
end

# `QuasiEwaldRbeInit`, formerly here, took an `ExTinyMD.SimulationInfo` and was
# removed while decoupling this package from ExTinyMD (see the phase-3 report).
# It is not a loss: it referenced `L` and `boundary`, neither ever defined in
# its own scope or passed as an argument, so calling it always raised
# `UndefVarError` -- it had no test coverage and was, as far as can be told,
# never callable in the form it existed in.
