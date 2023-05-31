@testset "compare the short range QEM with ICM" begin
    n_atoms = 100
    L = 100.0
    boundary = Q2dBoudary(L, L, 10.0)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms
        push!(atoms, Atom(mass = 1.0, charge = (-1.0)^i))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, 10.0), boundary; min_r = 1.0, temp = 1.0)

    interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
    loggers = [TempartureLogger(100)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    γ_1 = 2.0 * rand() - 1.0
    γ_2 = 2.0 * rand() - 1.0
    ϵ_0 = 1.0
    accuracy = 1e-4
    α = 0.5
    n_t = 30

    N_real = 0
    N_img = 100
    ICM_sys = IcmSys((γ_2, γ_1), (L, L, 10.0), N_real, N_img)
    ref_pos, ref_charge = IcmSysInit(ICM_sys, info.coords, [atom.charge for atom in atoms])
    force_icm = IcmForce(ICM_sys, info.coords, [atom.charge for atom in atoms], ref_pos, ref_charge)
    # force_self = IcmForce_self(ICM_sys, info.coords, [atom.charge for atom in atoms], ref_pos, ref_charge)

    gauss_para = GaussParameter(n_t)
    force_qem = [Point(0.0, 0.0, 0.0) for i in 1:n_atoms];
    for i in 1:n_atoms
        coord_1 = info.coords[i]
        q_1 = atoms[i].charge

        for j in i + 1:n_atoms    
            coord_2 = info.coords[j]
            q_2 = atoms[j].charge
            ge = GreensElement(γ_1, γ_2, coord_1[3], coord_2[3], sqrt((coord_1[1] - coord_2[1])^2 + (coord_1[2] - coord_2[2])^2), 10.0, α, 1e-4)
            force_i, force_j = QuasiEwald_Fs_pair(q_1, q_2, ϵ_0, ge, coord_1, coord_2, gauss_para; single_mode = true)
            force_qem[i] += force_i
            force_qem[j] += force_j
        end

        ge = GreensElement(γ_1, γ_2, coord_1[3], 10.0, α, 1e-4)
        force_qem[i] += QuasiEwald_Fs_self(q_1, ϵ_0, ge, gauss_para; single_mode = true)
    end

    for i in 1:n_atoms
        @test isapprox(force_icm[i][1], force_qem[i][1], atol = 1e-6)
        @test isapprox(force_icm[i][2], force_qem[i][2], atol = 1e-6)
        @test isapprox(force_icm[i][3], force_qem[i][3], atol = 1e-6)
    end

end