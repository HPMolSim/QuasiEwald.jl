@testset "framework-free plan API (Task 3)" begin
    # This is the standalone core query API this phase exists to add:
    # QuasiEwaldShortPlan/QuasiEwaldLongPlan constructed and queried directly
    # from plain arrays, no ExTinyMD type touched. (This file's *suite* does
    # `using ExTinyMD` via runtests.jl for other tests and for the ICM cross
    # check below, so it is not itself the "never loaded" proof -- that is
    # test/standalone.jl, added in Task 5. This test instead checks that the
    # plan API gives the physically correct answer.)
    n_atoms = 20
    L = 10.0
    Lz = 10.0

    # Positions/charges built directly, deliberately avoiding
    # ExTinyMD.SimulationInfo (see the phase plan's baseline-capture warning:
    # SimulationInfo consumes rand() internally).
    rng = Random.MersenneTwister(20260917)
    poses = [SVector(L * rand(rng), L * rand(rng), 1.0 + 8.0 * rand(rng)) for _ in 1:n_atoms]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n_atoms]

    for (γ_1, γ_2) in [(0.0, 0.0), (0.4, 0.5), (0.4, -0.5), (-0.4, -0.5)]
        ϵ_0 = 1.0
        n_t = 100
        accuracy = 1e-4
        α = 10.0
        r_c = 4.5   # r_c = 4.5 < L/2 = 5.0
        k_c = sqrt(-4 * α * log(accuracy))

        sp = QuasiEwaldShortPlan(γ_1, γ_2, ϵ_0, (L, L, Lz), false, accuracy, α, n_atoms, r_c, n_t)
        lp = QuasiEwaldLongPlan(γ_1, γ_2, ϵ_0, (L, L, Lz), false, accuracy, α, n_atoms, k_c, 0)

        Es = QuasiEwald.energy(sp, poses, charges)
        El = QuasiEwald.energy(lp, poses, charges)
        Fs = QuasiEwald.force(sp, poses, charges)
        Fl = QuasiEwald.force(lp, poses, charges)

        N_real, N_img = 100, 20
        ICM_sys = IcmSys((γ_2, γ_1), (L, L, Lz), N_real, N_img)
        ref_pos, ref_charge = IcmSysInit(ICM_sys, poses, charges)
        energy_icm = IcmEnergy(ICM_sys, poses, charges, ref_pos, ref_charge)
        force_icm = IcmForce(ICM_sys, poses, charges, ref_pos, ref_charge) ./ ϵ_0

        @testset "energy for γ = ($γ_1, $γ_2)" begin
            @test isapprox(energy_icm, Es + El, atol = 1e-2)
        end

        @testset "force for γ = ($γ_1, $γ_2)" begin
            for i in 1:n_atoms
                total = Fs[i] .+ Fl[i]
                error_i = sqrt(sum(abs2, total .- force_icm[i]))
                @test error_i < 1e-2
            end
        end

        @testset "force!/force in-place vs. allocating agree, for γ = ($γ_1, $γ_2)" begin
            F = [SVector(0.0, 0.0, 0.0) for _ in 1:n_atoms]
            QuasiEwald.force!(F, sp, poses, charges)
            @test F == Fs
        end
    end

    # Never mutate the caller's arrays.
    @testset "queries do not mutate poses/charges" begin
        poses_copy = deepcopy(poses)
        charges_copy = deepcopy(charges)
        sp = QuasiEwaldShortPlan(0.4, 0.5, 1.0, (L, L, Lz), false, 1e-4, 10.0, n_atoms, 4.5, 30)
        lp = QuasiEwaldLongPlan(0.4, 0.5, 1.0, (L, L, Lz), false, 1e-4, 10.0, n_atoms, sqrt(-4 * 10.0 * log(1e-4)), 0)
        QuasiEwald.energy(sp, poses, charges)
        QuasiEwald.force(sp, poses, charges)
        QuasiEwald.energy(lp, poses, charges)
        QuasiEwald.force(lp, poses, charges)
        @test poses == poses_copy
        @test charges == charges_copy
    end

    # Finite-difference self-consistency: F = -dE/dr for the plan's own
    # energy/force, independent of the (buggy, see report) old MDSys
    # adapter. Small n for speed; this is the test that would catch a
    # mass-division-style sign or scaling fault in the new core.
    @testset "force matches finite-difference of energy (short + long)" begin
        n = 6
        rng2 = Random.MersenneTwister(4)
        p = [SVector(L * rand(rng2), L * rand(rng2), 1.0 + 8.0 * rand(rng2)) for _ in 1:n]
        c = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
        # accuracy = 1e-8 and n_t = 100, rather than the 1e-4 / n_t = 30 used
        # elsewhere in this file, and the reason is the whole point of this test.
        # The short-range z-force is the one component that needs BOTH knobs
        # tightened, and they are not independent:
        #
        #   * `accuracy` sets where the k-space integrals are truncated. The
        #     resulting error floor is ~60 * accuracy.
        #   * `n_t` sets the Gauss quadrature order over that interval.
        #
        # Tightening `accuracy` ALONE makes the z-force worse, not better,
        # because a tighter truncation widens the interval that a fixed-order
        # rule has to cover. Measured max |F_z - (-dE/dz)| on this exact 6-particle
        # configuration:
        #
        #       n_t \ acc    1e-4      1e-6      1e-8     1e-10
        #       30         5.7e-6    7.5e-6    2.0e-5    3.8e-5   <- diverges
        #       60         6.3e-6    6.5e-8    4.1e-9    1.5e-8
        #       100        6.3e-6    6.5e-8    8.6e-10   1.4e-10
        #       400        6.3e-6    6.5e-8    8.6e-10   1.4e-10  <- converged
        #
        # So n_t = 30 (this file's default elsewhere) is simply not a converged
        # quadrature for the z-derivative, and at n_t >= 100 the error tracks
        # `accuracy` as it should. There is no formula defect in Fsz_point_core /
        # Fsz_gauss_core / dz_Gamma_* -- an earlier draft of this test excluded the
        # short-range z-component on the theory that there was one; the sweep above
        # is what settled it. With both knobs set, no component needs excluding.
        γ_1, γ_2, ϵ_0, accuracy, α, r_c, n_t = 0.4, 0.5, 1.0, 1e-8, 10.0, 4.5, 100
        k_c = sqrt(-4 * α * log(accuracy))
        sp = QuasiEwaldShortPlan(γ_1, γ_2, ϵ_0, (L, L, Lz), false, accuracy, α, n, r_c, n_t)
        lp = QuasiEwaldLongPlan(γ_1, γ_2, ϵ_0, (L, L, Lz), false, accuracy, α, n, k_c, 0)

        h = 1e-6
        for plan in (sp, lp)
            F = QuasiEwald.force(plan, p, c)
            for i in 1:n, d in 1:3
                pp = copy(p); pm = copy(p)
                vp = collect(p[i]); vp[d] += h
                vm = collect(p[i]); vm[d] -= h
                pp[i] = SVector{3,Float64}(vp...)
                pm[i] = SVector{3,Float64}(vm...)
                Ep = QuasiEwald.energy(plan, pp, c)
                Em = QuasiEwald.energy(plan, pm, c)
                fd = -(Ep - Em) / (2h)
                # No component is excluded: every plan, every direction.
                @test isapprox(fd, F[i][d], atol = 1e-8, rtol = 1e-3)
            end
        end
    end
end
