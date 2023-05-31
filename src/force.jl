function ExTinyMD.update_acceleration!(interaction::QuasiEwaldShortInteraction{T, TI}, neighborfinder::CellListDirQ2D{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    atoms = sys.atoms
    boundary = sys.boundary
    QuasiEwald_Fs!(interaction, neighborfinder, atoms, boundary, info.coords, info.acceleration)
    return nothing
end

function ExTinyMD.update_acceleration!(interaction::QuasiEwaldLongInteraction{T, TI}, neighbor::CellListDirQ2D{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    return nothing
end