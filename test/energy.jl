@testset "compare the energy with ICM method" begin
    n_atoms = 20
    L = 10.0
    boundary = ExTinyMD.Q2dBoundary(L, L, 10.0)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms
        push!(atoms, Atom(mass = 1.0, charge = (-1.0)^i))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, 10.0), boundary; min_r = 1.0, temp = 1.0)

    interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
    loggers = [TempartureLogger(100, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    for (γ_1, γ_2) in [(0.0, 0.0), (0.4, 0.5), (0.4, -0.5), (-0.4, -0.5)]
        ϵ_0 = 1.0
        n_t = 100

        N_real = 100
        N_img = 20
        ICM_sys = IcmSys((γ_2, γ_1), (L, L, 10.0), N_real, N_img)
        ref_pos, ref_charge = IcmSysInit(ICM_sys, info.coords, [atom.charge for atom in atoms])
        energy_icm = IcmEnergy(ICM_sys, info.coords, [atom.charge for atom in atoms], ref_pos, ref_charge)

        accuracy = 1e-4
        α = 10.0
        # r_c = (α * accuracy)^(-1/3)
        r_c = 4.5
        k_c = sqrt(-4 * α * log(accuracy))

        sortz = SortingFinder(info.coords)
        cellq2d = CellListDirQ2D(info, r_c, boundary, 1)
        interaction_short = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, r_c, n_t)
        interaction_long = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, k_c, 0)

        Es = QuasiEwald_Es(interaction_short, cellq2d, sys, info)
        El = QuasiEwald_El(interaction_long, sortz, sys, info)

        @test isapprox(energy_icm, Es + El, atol = 1e-2)
    end
end

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
        ra = RingAngles(k_0, Lx, Ly, Lz, α, k_c)
        split_sum = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c, ra)
        sort_sum = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)
        direct_sum = direct_sum_total(q, coords, L, γ_1, γ_2, ϵ_0, α, k_c)

        @test split_sum ≈ sort_sum
        @test split_sum ≈ direct_sum
        @test sort_sum ≈ direct_sum
    end
end