using QuasiEwald
using Test
using ExTinyMD
using SpecialFunctions
using StaticArrays
using Random

@testset "QuasiEwald.jl" begin
    include("Icm.jl")
    include("energy_short.jl")
    include("force.jl")
    include("energy.jl")
    include("force_long.jl")
    include("energy_long.jl")
    include("simulate.jl")
    include("plan.jl")
end
