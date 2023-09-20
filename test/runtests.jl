using QuasiEwald
using Test
using BenchmarkTools
using ExTinyMD

@testset "QuasiEwald.jl" begin
    include("test_Icm.jl")
    include("test_integrate.jl")
    include("test_force.jl")
    include("test_energy.jl")
    include("test_simulate.jl")
end
