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
    script_contaminated = """
    using ExTinyMD, QuasiEwald, StaticArrays
    @assert !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD"))
    print("OK")
    """
    proc = run(pipeline(
        `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script_contaminated`;
        stdout = devnull, stderr = devnull), wait = false)
    wait(proc)
    @test !success(proc)   # the @assert must trip -- confirms the guard is load-bearing
end
