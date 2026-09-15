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

    # ------------------------------------------------------------------------
    # Regression test for the ρ = 0 NaN (finding F1).
    #
    # `QuasiEwald_Fs_pair` builds the in-plane force from the RADIAL magnitude
    # `Fsr` times the unit vector (dx, dy)/ρ. At ρ == 0 that is 0/0: `Fsr`
    # itself vanishes there, so the unguarded quotient is NaN rather than a
    # genuine singularity. One NaN then propagates through the whole force
    # array and the whole trajectory, silently.
    #
    # Before the decoupling, ExTinyMD's `position_checkQ2D` returned its
    # all-zero sentinel for such a pair and every caller's `iszero(ρ_sq)`
    # guard skipped it -- which masked this by dropping the pair's real
    # z-force too. The explicit `ρ_sq ≥ r_c^2` test that replaced the
    # sentinel keeps the pair, so the guard has to live in the kernel.
    #
    # This asserts VALUES, not just `!isnan`: the in-plane components must be
    # exactly zero (symmetry: dx and dy are identically zero) and the
    # z-component must equal the dx -> 0 limit of the same configuration.
    # ------------------------------------------------------------------------
    @testset "ρ = 0 gives finite, correct forces (F1 regression)" begin
        L3 = (10.0, 10.0, 10.0)
        sp0 = QuasiEwaldShortPlan(0.4, 0.5, 1.0, L3, false, 1e-4, 10.0, 2, 4.5, 30)
        # r_c = 4.5 < min(Lx, Ly) / 2 = 5.0
        c2 = [1.0, -1.0]

        # Reference values, measured on the fixed code and independently
        # confirmed below against the dx -> 0 limit of the same pair.
        Fz_ref = 0.0001689588779540122
        E_ref = -0.4461255053088635

        # All three ways ρ = 0 arises. The minimum image reduces every one of
        # them to the same (ρ = 0, z_i = 2, z_j = 5) GreensElement, so all
        # three must give the identical answer.
        cases = (
            # (a) an exact (x, y) column -- what any lattice/column
            #     initialisation of a confined slab produces.
            ("shared (x, y) column", [SVector(3.0, 4.0, 2.0), SVector(3.0, 4.0, 5.0)]),
            # (b) in-plane separation an exact multiple of Lx, so the minimum
            #     image wraps to exactly zero.
            ("Δx = Lx (wraps to 0)", [SVector(3.0, 4.0, 2.0), SVector(13.0, 4.0, 5.0)]),
            ("Δx = -Lx (wraps to 0)", [SVector(3.0, 4.0, 2.0), SVector(-7.0, 4.0, 5.0)]),
            # (c) the same for Ly.
            ("Δy = Ly (wraps to 0)", [SVector(3.0, 4.0, 2.0), SVector(3.0, 14.0, 5.0)]),
        )

        for (label, poses0) in cases
            @testset "$label" begin
                F0 = QuasiEwald.force(sp0, poses0, c2)
                E0 = QuasiEwald.energy(sp0, poses0, c2)

                # Finiteness first: this is what regressed (NaN, not Inf).
                @test all(isfinite, (F0[1]..., F0[2]...))
                @test isfinite(E0)

                # In-plane components are *exactly* zero, not merely small:
                # dx and dy are identically zero, so symmetry admits no other
                # answer, and the guard sets them to `zero(T)`.
                @test F0[1][1] === 0.0
                @test F0[1][2] === 0.0
                @test F0[2][1] === 0.0
                @test F0[2][2] === 0.0

                # The z-component is finite and continuous through ρ = 0 and
                # must be left alone by the guard -- dropping the pair (the
                # pre-branch behaviour) would give 0 here instead.
                @test F0[1][3] ≈ Fz_ref rtol = 1e-12
                @test F0[2][3] ≈ -0.00015341820374866757 rtol = 1e-12
                @test !iszero(F0[1][3])
                @test E0 ≈ E_ref rtol = 1e-12
            end
        end

        # The dx -> 0 limit, computed independently of the ρ = 0 branch: only
        # dx > 0 configurations, which never reach the guard at all.
        @testset "ρ = 0 force agrees with the dx -> 0 limit" begin
            Fx_over_dx = Float64[]
            Fz_limit = 0.0
            for dx in (1e-1, 1e-2, 1e-3, 1e-6, 1e-9)
                poses_dx = [SVector(3.0, 4.0, 2.0), SVector(3.0 + dx, 4.0, 5.0)]
                Fdx = QuasiEwald.force(sp0, poses_dx, c2)
                push!(Fx_over_dx, Fdx[1][1] / dx)
                Fz_limit = Fdx[1][3]
            end

            # F_x vanishes LINEARLY in dx: F_x/dx is constant to 5 digits over
            # eight decades (measured 1.082565e-4, 1.088292e-4, 1.088349e-4,
            # 1.088350e-4, 1.088350e-4 at dx = 1e-1, 1e-2, 1e-3, 1e-6, 1e-9).
            # So the limit is 0, which is what the guard returns.
            @test Fx_over_dx[end] ≈ Fx_over_dx[end - 1] rtol = 1e-6
            @test Fx_over_dx[end] ≈ 1.088350e-4 rtol = 1e-5
            # and the raw F_x itself really does go to zero:
            @test abs(Fx_over_dx[end] * 1e-9) < 1e-12

            # F_z is smooth and even in dx, so F_z(dx) - F_z(0) = O(dx^2);
            # at dx = 1e-9 that correction is ~1e-18 relative, far below the
            # tolerance below (measured: bitwise equal in Float64).
            @test Fz_limit ≈ Fz_ref rtol = 1e-12

            poses_zero = [SVector(3.0, 4.0, 2.0), SVector(3.0, 4.0, 5.0)]
            @test QuasiEwald.force(sp0, poses_zero, c2)[1][3] ≈ Fz_limit rtol = 1e-12
        end

        # The original F1 reproduction: a third, well-separated particle must
        # not be contaminated -- a single NaN in the pair would propagate.
        @testset "a coincident pair does not poison the rest of the array" begin
            sp3 = QuasiEwaldShortPlan(0.4, 0.5, 1.0, L3, false, 1e-4, 10.0, 3, 4.5, 30)
            poses3 = [SVector(2.0, 3.0, 2.0), SVector(2.0, 3.0, 5.0), SVector(7.0, 8.0, 4.0)]
            c3 = [1.0, -1.0, 1.0]
            F3 = QuasiEwald.force(sp3, poses3, c3)
            @test all(isfinite, (F3[1]..., F3[2]..., F3[3]...))
            @test isfinite(QuasiEwald.energy(sp3, poses3, c3))
            @test F3[1][1] === 0.0 && F3[1][2] === 0.0
            @test F3[2][1] === 0.0 && F3[2][2] === 0.0
        end
    end
end
