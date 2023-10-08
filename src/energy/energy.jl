function ExTinyMD.energy(interaction::QuasiEwaldShortInteraction{T, TI}, neighborfinder::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer, T_NEIGHBOR<:ExTinyMD.AbstractNeighborFinder}
    update_finder!(neighborfinder, info)
    return QuasiEwald_Es(interaction, neighborfinder, sys, info)
end

function ExTinyMD.energy(interaction::QuasiEwaldLongInteraction{T, TI}, neighborfinder::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer, T_NEIGHBOR<:ExTinyMD.AbstractNeighborFinder}
    update_finder!(neighborfinder, info)
    return QuasiEwald_El(interaction, neighborfinder, sys, info)
end