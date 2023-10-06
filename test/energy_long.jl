@testset "energy_sum_total with direct_sum_total" begin
    n_atoms = 100
    q = 2 .* randn(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms

    Lx, Ly, Lz = (100.0, 100.0, 10.0)
    L = (Lx, Ly, Lz)

    boundary = ExTinyMD.Q2dBoundary(Lx, Ly, Lz)
    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms
        push!(atoms, Atom(mass = 1.0, charge = q[i]))
    end
    info = SimulationInfo(n_atoms, atoms, (0.0, Lx, 0.0, Ly, 0.0, Lz), boundary; min_r = 1.0, temp = 1.0)

    coords = info.coords
    z_coords = [coords[i][3] for i=1:n_atoms]
    z_list = sortperm(z_coords)
    for (γ_1, γ_2) in [(0.0, 0.0), (0.4, 0.5), (0.4, -0.5), (-0.4, -0.5)]
        ϵ_0 = 1.0
        k_c = 1.0
        α = 0.1
        sort_sum = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)
        direct_sum = direct_sum_total(q, coords, L, γ_1, γ_2, ϵ_0, α, k_c)

        @test sort_sum ≈ direct_sum
    end

    for (γ_1, γ_2) in [(10.0, 10.0), (-10.0, -10.0)]
        ϵ_0 = 1.0
        k_c = 1.0
        α = 0.1
        k_0 = log(γ_1 * γ_2) / (2 * Lz)
        ra = RingAngles(k_0, Lx, Ly, Lz, α, k_c, π/Lx)
        split_sum = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c, ra)
        sort_sum = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)
        direct_sum = direct_sum_total(q, coords, L, γ_1, γ_2, ϵ_0, α, k_c)

        @test split_sum ≈ sort_sum
        @test split_sum ≈ direct_sum
        @test sort_sum ≈ direct_sum
    end
end