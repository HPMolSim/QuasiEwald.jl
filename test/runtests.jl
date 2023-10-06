using QuasiEwald
using Test
using ExTinyMD
using SpecialFunctions

@testset "QuasiEwald.jl" begin
    include("Icm.jl")
    include("energy_short.jl")
    include("force.jl")
    include("energy.jl")
    include("force_long.jl")
    include("energy_long.jl")
    include("simulate.jl")
end
