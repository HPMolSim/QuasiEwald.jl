using QuasiEwald
using Test
using BenchmarkTools
using ExTinyMD

@testset "QuasiEwald.jl" begin
    include("test_Icm.jl")
end
