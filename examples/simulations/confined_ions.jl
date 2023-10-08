using ExTinyMD, QuasiEwald

begin
    n_atoms = 436
    n_atoms = Int64(round(n_atoms))
    L_x = 100.0
    L_y = 100.0
    L_z = 50.0
    L = (L_x, L_y, L_z)
    boundary = Q2dBoundary(L_x, L_y, L_z)
    atoms = Vector{Atom{Float64}}()

    for i in 1:218
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in 219:436
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    (γ_1, γ_2) = (0.95, -0.95)

    info = SimulationInfo(n_atoms, atoms, (0.0, L_x, 0.0, L_y, 0.5, L_z - 0.5), boundary; min_r = 2.0, temp = 1.0)

    ϵ_0 = 1.0

    accuracy = 1e-4
    α = 1.0
    k_c = sqrt(- 4 * α * log(accuracy))
    r_c = (α * accuracy)^(-1/3) / 2
    n_t = 30
    rbe_p = 50

    intershort = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, r_c, n_t)
    short_finder = CellListQ2D(info, r_c + 1.0, boundary, 100)
    interlong = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, k_c, rbe_p)
    long_finder = SortingFinder(info)

    interactions = [
        (LennardJones(), CellList3D(info, 4.5, boundary, 100)),
        (SubLennardJones(0.0, L_z; cutoff = 0.5, σ = 0.5), SubNeighborFinder(1.0, info, 0.0, L_z)), 
        (intershort, short_finder),
        (interlong, long_finder)
        ]

    loggers = [TempartureLogger(100, output = true), TrajectionLogger(step = 100, output = true)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    simulate!(simulator, sys, info, 100000)
end