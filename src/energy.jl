function ExTinyMD.energy(interaction::QuasiEwaldShortInteraction{T, TI}, neighbor::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer, T_NEIGHBOR<:ExTinyMD.AbstractNeighborFinder}
    return QuasiEwald_Es(interaction, neighbor.neighbor_list, sys.atoms, info.coords)
end

function ExTinyMD.energy(interaction::QuasiEwaldLongInteraction{T, TI}, neighbor::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer, T_NEIGHBOR<:ExTinyMD.AbstractNeighborFinder}
    return QuasiEwald_El(interaction, neighbor, sys, info)
end