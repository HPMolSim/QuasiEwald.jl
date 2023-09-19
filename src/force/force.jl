function ExTinyMD.update_acceleration!(interaction::QuasiEwaldShortInteraction{T, TI}, neighborfinder::CellListDirQ2D{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    update_finder!(neighborfinder, info)
    QuasiEwald_Fs!(interaction, neighborfinder, sys.atoms, sys.boundary, info.coords, info.acceleration)
    return nothing
end

function ExTinyMD.update_acceleration!(interaction::QuasiEwaldLongInteraction{T, TI}, neighborfinder::SortingFinder{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    update_finder!(neighborfinder, info)
    QuasiEwald_Fl!(interaction, neighborfinder, sys.atoms, sys.boundary, info.coords, info.acceleration)
    return nothing
end