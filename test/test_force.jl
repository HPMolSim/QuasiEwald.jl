@testset "compare the force with ICM method" begin
    n_atoms = 20
    L = 10.0
    boundary = ExTinyMD.Q2dBoundary(L, L, 10.0)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms
        push!(atoms, Atom(mass = 1.0, charge = (-1.0)^i))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, 10.0), boundary; min_r = 1.0, temp = 1.0)

    interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
    loggers = []
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
        force_icm = IcmForce(ICM_sys, info.coords, [atom.charge for atom in atoms], ref_pos, ref_charge) ./ ϵ_0;

        accuracy = 1e-4
        α = 1.0
        # r_c = (α * accuracy)^(-1/3)
        r_c = 4.5
        k_c = sqrt(-4 * α * log(accuracy))

        sortz = SortingFinder(info.coords)
        cellq2d = CellListDirQ2D(info, r_c, boundary, 1)
        interaction_short = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, r_c, n_t)
        interaction_long = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, k_c, 0)

        force_qem = [Point(0.0, 0.0, 0.0) for i in 1:n_atoms]
        QuasiEwald_Fs!(interaction_short, cellq2d, atoms, boundary, info.coords, force_qem)
        force_qem_long = [Point(0.0, 0.0, 0.0) for i in 1:n_atoms]
        QuasiEwald_Fl!(interaction_long, sortz, atoms, boundary, info.coords, force_qem)

        for i in 1:n_atoms
            error_i = force_icm[i] - force_qem[i]
            for j in 1:3
                @test abs(error_i[j]) < 1e-3
            end
        end
    end
end