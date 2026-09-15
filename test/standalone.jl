@testset "core works without ExTinyMD" begin
    # The requirement this whole phase exists for: QuasiEwald's core query API
    # must work with ExTinyMD never loaded. A @testset inside the normal suite
    # does not prove that on its own -- test/runtests.jl itself does `using
    # ExTinyMD` (for the adapter tests), so by the time this testset runs,
    # ExTinyMD is already loaded in *this* process. The only way to prove the
    # core doesn't need it is to run in a fresh process that never imports it.
    # Modelled on ParticleMeshEwald's test/standalone.jl.
    script = """
    using QuasiEwald, StaticArrays
    @assert !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD"))

    n = 24
    L = (12.0, 12.0, 8.0)
    poses = [SVector(rand()*L[1], rand()*L[2], 1.0 + rand()*(L[3]-2.0)) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    ϵ_0 = 1.0
    accuracy = 1e-4
    α = 1.0
    r_c = 4.5      # r_c = 4.5 < min(Lx,Ly)/2 = 6.0
    k_c = sqrt(-4 * α * log(accuracy))
    n_t = 30

    sp = QuasiEwaldShortPlan(0.4, 0.5, ϵ_0, L, false, accuracy, α, n, r_c, n_t)
    lp = QuasiEwaldLongPlan(0.4, 0.5, ϵ_0, L, false, accuracy, α, n, k_c, 0)

    Es = QuasiEwald.energy(sp, poses, charges)
    El = QuasiEwald.energy(lp, poses, charges)
    Fs = QuasiEwald.force(sp, poses, charges)
    Fl = QuasiEwald.force(lp, poses, charges)
    @assert isfinite(Es) && isfinite(El)
    @assert all(isfinite, Fs[i][d] for i in 1:n, d in 1:3)
    @assert all(isfinite, Fl[i][d] for i in 1:n, d in 1:3)

    # the ICM path is framework-free too
    sys = IcmSys((0.4, 0.5), L, 20, 8)
    ref_pos, ref_charge = IcmSysInit(sys, poses, charges)
    E_icm = IcmEnergy(sys, poses, charges, ref_pos, ref_charge)
    F_icm = IcmForce(sys, poses, charges, ref_pos, ref_charge)
    @assert isfinite(E_icm)
    @assert all(isfinite, F_icm[i][d] for i in 1:n, d in 1:3)

    print("OK")
    """
    out = read(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script`, String)
    @test out == "OK"

    # Verify the check above is not vacuous: it must actually fail if
    # ExTinyMD is loaded first. Across this phase's predecessors, eleven
    # specifications turned out to be tests that passed for the wrong
    # reason; this is the one most at risk of that, so it is asserted
    # against directly rather than assumed.
    #
    # `@test !success(proc)` alone is NOT enough, and that is the whole point:
    # it is true for any non-zero exit -- ExTinyMD missing from the test
    # environment, an unsatisfiable resolve, a failed precompile, a typo in
    # the heredoc. Every one of those would make this check pass while proving
    # nothing about whether the `@assert` on `Base.loaded_modules` is
    # load-bearing. So the subprocess's streams are captured and the failure
    # is pinned to the intended cause: a marker printed after `using` and
    # before the assert must appear on stdout (so the loads all succeeded),
    # the marker after the assert must NOT appear (so the assert is what
    # stopped it), and stderr must carry the AssertionError naming the guard's
    # own expression.
    script_contaminated = """
    using ExTinyMD, QuasiEwald, StaticArrays
    print("LOADS_OK;")
    @assert !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD"))
    print("ASSERT_DID_NOT_TRIP;")
    """
    outfile, errfile = tempname(), tempname()
    proc = run(pipeline(
        `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script_contaminated`;
        stdout = outfile, stderr = errfile), wait = false)
    wait(proc)
    out_text = read(outfile, String)
    err_text = read(errfile, String)
    rm(outfile, force = true)
    rm(errfile, force = true)

    @test !success(proc)
    # `using ExTinyMD, QuasiEwald, StaticArrays` really did succeed, so the
    # non-zero exit is not a missing package, a resolve failure or a
    # precompilation error.
    @test occursin("LOADS_OK;", out_text)
    # and execution stopped at the assert, not after it.
    @test !occursin("ASSERT_DID_NOT_TRIP;", out_text)
    # and it stopped for exactly the intended reason.
    @test occursin("AssertionError", err_text)
    @test occursin("loaded_modules", err_text)

    # The dispatcher's two error cases (finding F6b). `Base.get_extension`
    # returns `nothing` both when ExTinyMD was never loaded and when it IS
    # loaded but the extension failed to precompile; reporting the second as
    # the first tells the user to run a `using` they have already run, while
    # the real error has scrolled past. Only the first case is reachable from
    # a test (the second needs a deliberately broken extension), so it is the
    # one asserted -- with the "loaded but the extension failed" wording
    # explicitly excluded, which is what pins the branch.
    script_no_extinymd = """
    using QuasiEwald
    try
        QuasiEwaldShortInteraction(0.4, 0.5, 1.0, (10.0, 10.0, 10.0), false, 1e-4, 1.0, 2, 4.5, 30)
        print("NO_ERROR")
    catch e
        print(sprint(showerror, e))
    end
    """
    msg = read(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script_no_extinymd`, String)
    @test occursin("requires ExTinyMD to be loaded", msg)
    @test occursin("QuasiEwaldExTinyMDExt", msg)
    @test !occursin("retry_load_extensions", msg)   # that is the *other* branch
end
