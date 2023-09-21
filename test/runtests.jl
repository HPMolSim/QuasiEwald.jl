using QuasiEwald
using Test
using BenchmarkTools
using ExTinyMD
using SpecialFunctions

@testset "QuasiEwald.jl" begin
    include("test_Icm.jl")
    include("test_integral.jl")
    include("test_force.jl")
    include("test_energy.jl")
    # include("test_simulate.jl")
end
