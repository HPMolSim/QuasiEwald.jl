@testset "Test QuasiEwald with ExTinyMD" begin
    n_atoms = 100
    n_atoms = Int64(round(n_atoms))
    L_x = 100.0
    L_y = 100.0
    L_z = 10.0
    L = (L_x, L_y, L_z)
    boundary = Q2dBoundary(L_x, L_y, L_z)
    atoms = Vector{Atom{Float64}}()

    for i in 1:50
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in 51:100
        push!(atoms, Atom(type = 2, mass = 1.0, charge = (-1.0)))
    end

    for (γ_1, γ_2) in [(0.0, 0.0), (0.4, 0.5), (0.4, -0.5), (-0.4, -0.5), (10.0, 10.0), (10.0, -10.0)]
        @testset "testing MD run for γ = ($γ_1, $γ_2)" begin

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

            loggers = [TemperatureLogger(100, output = false), TrajectoryLogger(step = 100, output = false)]
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
end