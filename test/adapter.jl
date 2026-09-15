@testset "ExTinyMD adapter (Task 4): wrapper under load" begin
    # This is the test the phase plan calls out as mattering more than it
    # looks: ParticleMeshEwald never exercised the §4.3a wrapper pattern
    # (it has no forces, so its extension only ever supplies ExTinyMD.energy
    # called directly). QuasiEwald's wrappers ARE placed in sys.interactions
    # and driven through simulate! below -- the only thing that catches the
    # id/slot, mass-division and accumulation faults that a per-call test
    # cannot.

    function _charged_system(n, L, Lz)
        boundary = Q2dBoundary(L, L, Lz)
        atoms = Atom{Float64}[]
        for i in 1:(n ÷ 2)
            push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
        end
        for i in (n ÷ 2 + 1):n
            push!(atoms, Atom(type = 2, mass = 1.0, charge = -1.0))
        end
        info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 1.0, Lz - 1.0), boundary; min_r = 1.0, temp = 1.0)
        info.running_step = 1
        return boundary, atoms, info
    end

    @testset "adapter runs inside simulate! with bounded energy drift" begin
        Random.seed!(20260938)
        n, L, Lz = 20, 12.0, 10.0
        boundary, atoms, info = _charged_system(n, L, Lz)

        ϵ_0 = 1.0
        accuracy = 1e-4
        α = 1.0
        r_c = 4.5   # r_c = 4.5 < min(L,L)/2 = 6.0
        k_c = sqrt(-4 * α * log(accuracy))
        n_t = 30

        intershort = QuasiEwaldShortInteraction(0.0, 0.0, ϵ_0, (L, L, Lz), false, accuracy, α, n, r_c, n_t)
        short_finder = CellListQ2D(info, r_c, boundary, 1)
        interlong = QuasiEwaldLongInteraction(0.0, 0.0, ϵ_0, (L, L, Lz), false, accuracy, α, n, k_c, 0)
        long_finder = SortingFinder(info)

        sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                    interactions = [(intershort, short_finder), (interlong, long_finder)],
                    loggers = [TemperatureLogger(200; output = false)],
                    simulator = VerletProcess(dt = 1e-4))

        E0 = ExTinyMD.energy(intershort, short_finder, sys, info) +
             ExTinyMD.energy(interlong, long_finder, sys, info)
        simulate!(sys.simulator, sys, info, 200)
        E1 = ExTinyMD.energy(intershort, short_finder, sys, info) +
             ExTinyMD.energy(interlong, long_finder, sys, info)

        @test isfinite(E1)
        # Microcanonical Verlet at this dt should not let the electrostatic
        # energy run away; a sign error, a mass-division fault, or an
        # accumulation fault (double-counting or dropping pairs) shows up as
        # an unbounded value. This bound is generous (this is a correctness
        # smoke test, not a symplectic-integrator accuracy benchmark).
        @test abs(E1 - E0) < 0.1 * max(abs(E0), 1.0)
    end

    @testset "adapter is correct when slot order differs from id order" begin
        # In stock ExTinyMD `particle_info[i].id == i`, so a gather that
        # confuses slot with id looks correct forever. Permute the mapping so
        # the two differ, following ExTinyMD's own test_adapter.jl pattern:
        # give every id a distinct mass and charge (so a mix-up cannot
        # cancel out), reverse the slot order, and keep ids attached to
        # their particles via info.id_dict.
        Random.seed!(20260939)
        n, L, Lz = 12, 10.0, 8.0
        boundary, atoms, info = _charged_system(n, L, Lz)

        atoms = [Atom(type = a.type, mass = 1.0 + 0.1 * i,
                     charge = (isodd(i) ? 1.0 : -1.0) * (1 + 0.01 * i))
                for (i, a) in enumerate(atoms)]

        reverse!(info.particle_info)
        for i in eachindex(info.particle_info)
            info.id_dict[info.particle_info[i].id] = i
        end
        @test info.particle_info[1].id != 1   # the mapping really is permuted

        ϵ_0 = 1.0
        accuracy = 1e-4
        α = 1.0
        r_c = 3.5   # r_c = 3.5 < L/2 = 5.0
        k_c = sqrt(-4 * α * log(accuracy))
        n_t = 30

        intershort = QuasiEwaldShortInteraction(0.4, 0.5, ϵ_0, (L, L, Lz), false, accuracy, α, n, r_c, n_t)
        short_finder = CellListQ2D(info, r_c, boundary, 1)
        interlong = QuasiEwaldLongInteraction(0.4, 0.5, ϵ_0, (L, L, Lz), false, accuracy, α, n, k_c, 0)
        long_finder = SortingFinder(info)

        sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                    interactions = [(intershort, short_finder), (interlong, long_finder)],
                    loggers = [TemperatureLogger(100; output = false)],
                    simulator = VerletProcess(dt = 0.001))

        poses = [SVector(p.position[1], p.position[2], p.position[3]) for p in info.particle_info]
        charges = [atoms[p.id].charge for p in info.particle_info]

        Es_ref = QuasiEwald.energy(intershort.plan, poses, charges)
        El_ref = QuasiEwald.energy(interlong.plan, poses, charges)

        @test isapprox(ExTinyMD.energy(intershort, short_finder, sys, info), Es_ref, rtol = 1e-10)
        @test isapprox(ExTinyMD.energy(interlong, long_finder, sys, info), El_ref, rtol = 1e-10)

        Fs_ref = QuasiEwald.force(intershort.plan, poses, charges)
        Fl_ref = QuasiEwald.force(interlong.plan, poses, charges)

        for p in info.particle_info
            p.acceleration = Point(0.0, 0.0, 0.0)
        end
        ExTinyMD.update_acceleration!(intershort, short_finder, sys, info)
        ExTinyMD.update_acceleration!(interlong, long_finder, sys, info)

        for (slot, p) in enumerate(info.particle_info)
            m = atoms[p.id].mass
            for d in 1:3
                @test isapprox(p.acceleration[d], (Fs_ref[slot][d] + Fl_ref[slot][d]) / m, rtol = 1e-10)
            end
        end
    end
end
