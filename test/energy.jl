@testset "compare the energy with ICM method" begin
    n_atoms = 20
    L = 10.0
    boundary = ExTinyMD.Q2dBoundary(L, L, 10.0)

    atoms = Vector{Atom{Float64}}()
    for i in 1:10
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in 11:20
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
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

        accuracy = 1e-4
        α = 10.0
        # r_c = (α * accuracy)^(-1/3)
        r_c = 4.5
        k_c = sqrt(-4 * α * log(accuracy))

        sortz = SortingFinder(info)
        cellq2d = CellListDirQ2D(info, r_c, boundary, 1)
        interaction_short = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, r_c, n_t)
        interaction_long = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, k_c, 0)

        Es = ExTinyMD.energy(interaction_short, cellq2d, sys, info)
        El = ExTinyMD.energy(interaction_long, sortz, sys, info)

        N_real = 100
        N_img = 20
        ICM_sys = IcmSys((γ_2, γ_1), (L, L, 10.0), N_real, N_img)
        coords = [p_info.position for p_info in info.particle_info]
        charge = [atoms[p_info.id].charge for p_info in info.particle_info]
        ref_pos, ref_charge = IcmSysInit(ICM_sys, coords, charge)
        energy_icm = IcmEnergy(ICM_sys, coords, charge, ref_pos, ref_charge)

        @test isapprox(energy_icm, Es + El, atol = 1e-2)
    end
end