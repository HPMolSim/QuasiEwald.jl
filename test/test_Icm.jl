@testset "Image charge method work test" begin
    N = 10
    γ_up = 0.10
    γ_down = - 0.15

    L = (1.0, 1.0, 1.0)
    position = [Point(L[1] * rand(), L[2] * rand(), L[3] * rand()) for _=1:N]
    charge = [(-1.0)^i for i in 1:N]


    N_img = 5
    N_real = 5

    sys = IcmSys((γ_up, γ_down), L, N_real, N_img)
    ref_pos, ref_charge = IcmSysInit(sys, position, charge)

    energy_jl = IcmEnergy(sys, position, charge, ref_pos, ref_charge)
    force_jl = IcmForce(sys, position, charge, ref_pos, ref_charge)
    @test eltype(force_jl) <: Point{3, Float64}
end