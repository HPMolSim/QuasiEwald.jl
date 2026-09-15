@testset "ExTinyMD adapter (Task 4): wrapper under load" begin
    # This is the test the phase plan calls out as mattering more than it
    # looks: ParticleMeshEwald never exercised the §4.3a wrapper pattern
    # (it has no forces, so its extension only ever supplies ExTinyMD.energy
    # called directly). QuasiEwald's wrappers ARE placed in sys.interactions
    # and driven through simulate! below -- the only thing that catches the
    # id/slot, mass-division and accumulation faults that a per-call test
    # cannot.

    function _charged_system(n, L, Lz; temp = 1.0)
        boundary = Q2dBoundary(L, L, Lz)
        atoms = Atom{Float64}[]
        for i in 1:(n ÷ 2)
            push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
        end
        for i in (n ÷ 2 + 1):n
            push!(atoms, Atom(type = 2, mass = 1.0, charge = -1.0))
        end
        info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 1.0, Lz - 1.0), boundary; min_r = 1.0, temp = temp)
        info.running_step = 1
        return boundary, atoms, info
    end

    # Give every id a distinct mass and a distinct charge magnitude, so that a
    # mass-division or id/slot fault cannot cancel out. Velocities were drawn
    # against the uniform-mass atoms above; that only changes the initial
    # condition, not the conservation law being tested.
    _nonuniform(atoms, n) = [Atom(type = a.type, mass = 0.5 + 0.15 * i,
                                  charge = (i <= n ÷ 2 ? 1.0 : -1.0) * (1 + 0.02 * i))
                             for (i, a) in enumerate(atoms)]

    _kinetic(atoms, info) = sum(0.5 * atoms[p.id].mass *
                                (p.velocity[1]^2 + p.velocity[2]^2 + p.velocity[3]^2)
                                for p in info.particle_info)

    @testset "adapter conserves total energy inside simulate! (F3)" begin
        # What this used to assert, and why it was replaced.
        #
        # The old assertion was `abs(E1 - E0) < 0.1 * max(abs(E0), 1.0)` on the
        # ELECTROSTATIC energy alone, with a measured drift of 7.6e-4 against a
        # bound of 1e-1 -- a 130x margin. Worse, that quantity cannot be made
        # to fail by tightening the bound, because it barely responds to a
        # force fault at all: with an exact 2x electrostatic force error (both
        # (interaction, finder) pairs listed twice in sys.interactions) the
        # electrostatic drift measured 9.41e-4 against the correct 9.46e-4 --
        # the same to two digits. E_elec is ~1% of the total energy here, so
        # the trajectory's electrostatic energy after N steps is set by thermal
        # motion, not by whether the force is right.
        #
        # The quantity that DOES respond is the conserved one. `VerletProcess`
        # with the default `NoThermoStat` is symplectic, and electrostatics is
        # the only interaction in `sys`, so KE + E_elec is conserved to
        # O(dt^2). Double the force and the integrator is conserving
        # KE + 2*E_elec instead, so KE + E_elec drifts by the electrostatic
        # work -- a first-order, not a second-order, effect.
        #
        # Measured |Δ(KE + E_elec)| over 500 steps at dt = 5e-3, γ = (0.4, 0.5),
        # non-uniform mass, temp = 0.05:
        #
        #   correct forces      1.518e-3   (worst 1.52e-3 over 8 different seeds)
        #   2x force fault      1.115e-1   (smallest 8.3e-3 over the same 8 seeds)
        #
        # Bound set to 1e-2: 6.6x above the measured correct drift (and above
        # the worst of the eight seeds), 11x below the faulted drift at this
        # test's own seed. Verified to fail on the 2x fault and pass without it.
        #
        # γ = (0.4, 0.5) rather than the old (0.0, 0.0): a zero γ pair zeroes
        # `force_k_sum_2!` (scales with γ_1*γ_2), `force_k_sum_3!` (γ_1) and
        # `force_k_sum_4!` (γ_2) outright, leaving only `force_k_sum_1!`
        # exercised. temp = 0.05 rather than 1.0 keeps the slab from flying
        # apart over 500 steps at this dt, with no thermostat to hold it.
        Random.seed!(20260938)
        n, L, Lz = 20, 12.0, 10.0
        boundary, atoms, info = _charged_system(n, L, Lz; temp = 0.05)
        atoms = _nonuniform(atoms, n)

        ϵ_0 = 1.0
        accuracy = 1e-4
        α = 1.0
        r_c = 4.5   # r_c = 4.5 < min(L,L)/2 = 6.0
        k_c = sqrt(-4 * α * log(accuracy))
        n_t = 30

        intershort = QuasiEwaldShortInteraction(0.4, 0.5, ϵ_0, (L, L, Lz), false, accuracy, α, n, r_c, n_t)
        short_finder = CellListQ2D(info, r_c, boundary, 1)
        interlong = QuasiEwaldLongInteraction(0.4, 0.5, ϵ_0, (L, L, Lz), false, accuracy, α, n, k_c, 0)
        long_finder = SortingFinder(info)

        sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                    interactions = [(intershort, short_finder), (interlong, long_finder)],
                    loggers = [TemperatureLogger(10^9; output = false)],
                    simulator = VerletProcess(dt = 5e-3))

        Eel0 = ExTinyMD.energy(intershort, short_finder, sys, info) +
               ExTinyMD.energy(interlong, long_finder, sys, info)
        K0 = _kinetic(atoms, info)
        simulate!(sys.simulator, sys, info, 500)
        Eel1 = ExTinyMD.energy(intershort, short_finder, sys, info) +
               ExTinyMD.energy(interlong, long_finder, sys, info)
        K1 = _kinetic(atoms, info)

        @test isfinite(Eel1)
        @test isfinite(K1)
        # The load-bearing assertion: the forces the integrator used must be
        # the gradient of the energy these same wrappers report. A sign error,
        # a scaling error, a dropped or double-counted mass division, or a
        # dropped/double-counted pair all break that equality and show up here.
        @test abs((K1 + Eel1) - (K0 + Eel0)) < 1e-2
        # The electrostatic energy itself must at least stay in range -- this
        # one is a smoke check only, for the reason given above.
        @test abs(Eel1 - Eel0) < 1.0
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

    @testset "short-range wrapper rejects 3-D neighbour finders (F7)" begin
        # `_finder_list(f) = f.neighbor_list` with an untyped `neighborfinder`
        # accepted CellList3D/CellListDir3D, which build their lists from 3-D
        # distances and therefore drop pairs that are close in plane but far
        # apart in z -- so the quasi-2D short-range sum came back silently too
        # small. Pre-decoupling those were a MethodError; the fallback method
        # restores a loud failure.
        Random.seed!(20260940)
        n, L, Lz = 12, 10.0, 8.0
        boundary, atoms, info = _charged_system(n, L, Lz)

        ϵ_0, accuracy, α, n_t = 1.0, 1e-4, 1.0, 30
        r_c = 3.5   # r_c = 3.5 < min(L,L)/2 = 5.0

        intershort = QuasiEwaldShortInteraction(0.4, 0.5, ϵ_0, (L, L, Lz), false, accuracy, α, n, r_c, n_t)
        sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                    interactions = [(intershort, CellListQ2D(info, r_c, boundary, 1))],
                    loggers = [TemperatureLogger(100; output = false)],
                    simulator = VerletProcess(dt = 0.001))

        for bad in (CellList3D(info, r_c, boundary, 1), CellListDir3D(info, r_c, boundary, 1))
            @test_throws ArgumentError ExTinyMD.energy(intershort, bad, sys, info)
            @test_throws ArgumentError ExTinyMD.update_acceleration!(intershort, bad, sys, info)
        end

        # The supported finders all still work, and -- since the plan always
        # recomputes the true in-plane distance from the candidate list -- the
        # two quasi-2D finders and the O(n^2) fallback must agree exactly.
        E_q2d = ExTinyMD.energy(intershort, CellListQ2D(info, r_c, boundary, 1), sys, info)
        E_dirq2d = ExTinyMD.energy(intershort, CellListDirQ2D(info, r_c, boundary, 1), sys, info)
        E_none = ExTinyMD.energy(intershort, NoNeighborFinder(), sys, info)
        @test isapprox(E_q2d, E_none, rtol = 1e-12)
        @test isapprox(E_dirq2d, E_none, rtol = 1e-12)

        # And the message has to say what to use, not just that it failed.
        msg = try
            ExTinyMD.energy(intershort, CellList3D(info, r_c, boundary, 1), sys, info)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("CellListQ2D", msg)
        @test occursin("CellList3D", msg)
    end
    @testset "the three preserved names are functions, not types (F6a)" begin
        # Documented breaking change: `QuasiEwaldShortInteraction` and friends
        # used to be structs and are now dispatcher functions, so construction
        # is unchanged but every type-position use of the bare name raises a
        # TypeError. This test is the executable form of that note in the
        # README and in the dispatcher's docstring -- if the pattern is ever
        # changed back, this is what says the documentation has gone stale.
        Random.seed!(20260941)
        n, L, Lz = 6, 10.0, 8.0
        boundary, atoms, info = _charged_system(n, L, Lz)
        α, accuracy = 1.0, 1e-4
        inter = QuasiEwaldShortInteraction(0.4, 0.5, 1.0, (L, L, Lz), false, accuracy, α, n, 3.5, 30)
        finder = SortingFinder(info)

        # Construction works exactly as before, and the result still subtypes
        # ExTinyMD's abstract types.
        @test inter isa ExTinyMD.AbstractInteraction
        @test finder isa ExTinyMD.AbstractNeighborFinder

        # But the exported names are not types.
        @test !(QuasiEwaldShortInteraction isa Type)
        @test !(QuasiEwaldLongInteraction isa Type)
        @test !(SortingFinder isa Type)
        @test_throws TypeError (inter isa QuasiEwaldShortInteraction)
        @test_throws TypeError (finder isa SortingFinder)

        # The documented workaround.
        ext = Base.get_extension(QuasiEwald, :QuasiEwaldExTinyMDExt)
        @test ext !== nothing
        @test inter isa ext.QuasiEwaldShortInteraction
        @test finder isa ext.SortingFinder
        @test QuasiEwaldLongInteraction(0.4, 0.5, 1.0, (L, L, Lz), false, accuracy, α, n,
                                        sqrt(-4 * α * log(accuracy)), 0) isa ext.QuasiEwaldLongInteraction
    end
end
