using Plots, ExTinyMD, QuasiEwald

begin
    n_atoms = 2
    L = 180.0
    boundary = ExTinyMD.Q2dBoundary(L, L, 10.0)

    atoms = [Atom(type = 1, mass = 1.0, charge = 1.0), Atom(type = 2, mass = 1.0, charge = -1.0)]

    sys = MDSys(
                n_atoms = n_atoms,
                atoms = atoms,
                boundary = boundary,
                interactions = [(NoInteraction(), NoNeighborFinder(n_atoms))],
                loggers = [TrajectionLogger(step = 100, output = false)],
                simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
            )

    Force_x = Vector{Vector{Float64}}()
    coord_1 = Point(50.0, 50.0, 1.0)
    X = 0.1:0.1:40.0

    for (γ_1, γ_2) in [(0.0, 0.0), (0.95, 0.95), (-0.95, -0.95), (10.0, 10.0), (-10.0, -10.0)]
        ϵ_0 = 1.0
        n_t = 30

        accuracy = 1e-4
        α = 1.0
        r_c = 15.0
        k_c = sqrt(-4 * α * log(accuracy))

        force_x = Vector{Float64}()

        for x_2 in X
            info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, 10.0), boundary; min_r = 1.0, temp = 1.0)
            coord_2 = Point(50.0 + x_2, 50.0, 1.01)

            info.particle_info[1].position = coord_1
            info.particle_info[2].position = coord_2

            sortz = SortingFinder(info)
            cellq2d = CellListQ2D(info, r_c, boundary, 1)
            interaction_short = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, r_c, n_t)
            interaction_long = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, k_c, 0)
    
            force_qem = [Point(0.0, 0.0, 0.0) for i in 1:n_atoms]
            QuasiEwald_Fs!(interaction_short, cellq2d, sys, info)
            QuasiEwald_Fl!(interaction_long, sortz, sys, info)
            push!(force_x, info.particle_info[1].acceleration[1])
        end
        push!(Force_x, force_x)
    end

    plot(dpi = 300, size = (800, 600), legend = :topright, xlabel = "x", ylabel = "force_x")
    for (γ, force_x) in zip([0.0, 0.95, -0.95, 10.0, -10.0], Force_x)
        plot!(X, force_x, label = "γ = " * string(γ), ylim = [-0.06, 0.06])
    end
    savefig("force_x.png")
end