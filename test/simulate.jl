@testset "Test QuasiEwald with ExTinyMD" begin
    n_atoms = 100
    n_atoms = Int64(round(n_atoms))
    L_x = 100.0
    L_y = 100.0
    L_z = 30.0
    L = (L_x, L_y, L_z)
    boundary = Q2dBoundary(L_x, L_y, L_z)
    atoms = Vector{Atom{Float64}}()

    for i in 1:n_atoms
        push!(atoms, Atom(mass = 1.0, charge = (-1.0)^i))
    end

    for (γ_1, γ_2) in [(0.0, 0.0), (0.4, 0.5), (0.4, -0.5), (-0.4, -0.5), (10.0, 10.0), (10.0, -10.0)]

        info = SimulationInfo(n_atoms, atoms, (0.0, L_x, 0.0, L_y, 0.5, L_z - 0.5), boundary; min_r = 2.0, temp = 1.0)
        coords = info.coords;
        z_coords = [coord[3] for coord in coords];
        z_list = sortperm(z_coords);

        ϵ_0 = 1.0

        accuracy = 1e-4
        α = 1.0
        k_c = sqrt(- 4 * α * log(accuracy))
        r_c = (α * accuracy)^(-1/3) / 2
        n_t = 30
        rbe_p = 50

        intershort = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, r_c, n_t)
        short_finder = CellListDirQ2D(info, r_c + 1.0, boundary, 100)
        interlong = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, k_c, rbe_p)
        long_finder = SortingFinder(coords)

        interactions = [
            (LennardJones(), CellListDir3D(info, 4.5, boundary, 100)),
            (SubLennardJones(0.0, L_z; cutoff = 0.5, σ = 0.5), SubNeighborFinder(1.0, info.coords, 0.0, L_z)), 
            (intershort, short_finder),
            (interlong, long_finder)
            ]

        loggers = [TempartureLogger(100, output = false), TrajectionLogger(info, 100, output = false)]
        simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

        sys = MDSys(
            n_atoms = n_atoms,
            atoms = atoms,
            boundary = boundary,
            interactions = interactions,
            loggers = loggers,
            simulator = simulator
        )
        simulate!(simulator, sys, info, 200)

        @test info.running_step == 200
    end
end