using QuasiEwald
using Test
using ExTinyMD
using SpecialFunctions

@testset "QuasiEwald.jl" begin
    include("Icm.jl")
    include("integral.jl")
    include("force.jl")
    include("energy.jl")
    # include("test_simulate.jl")
end
