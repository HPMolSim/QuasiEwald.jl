n_atoms = 1000
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

γ_1 = 0.0
γ_2 = - 0.0
ϵ_0 = 1.0
accuracy = 1e-4
α = 0.5
n_t = 100

Es = []
rc_array = 1.0:.5:10.0
for r_c in rc_array
    cellq2d = CellListDirQ2D(info, r_c, boundary, 1)
    interaction_short = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, r_c, n_t)
    @show r_c, length(cellq2d.neighbor_list)
    @time push!(Es, QuasiEwald_Es(interaction_short, cellq2d, sys, info))
end
plot(rc_array, Es)

El = []
kc_array = 1.0:.4:5.0
sortz = SortingFinder(info.coords)
for k_c in kc_array
    interaction_long = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, k_c, 0)
    @show k_c
    @time push!(El, QuasiEwald_El(interaction_long, sortz, sys, info))
end
plot(kc_array, El)

El[end] + Es[end]

N_real = 0
N_img = 100
ICM_sys = IcmSys((γ_2, γ_1), (L, L, 10.0), N_real, N_img)
ref_pos, ref_charge = IcmSysInit(ICM_sys, info.coords, [atom.charge for atom in atoms])
energy_jl = IcmEnergy(ICM_sys, info.coords, [atom.charge for atom in atoms], ref_pos, ref_charge)

accuracy = 1e-2
α = (π^2/(L^4 * log10(1/accuracy) * accuracy^(2/3)) * 18 * n_t / 16 * n_atoms)^0.6
r_c = (α * accuracy)^(-1/3)
k_c = sqrt(-4 * α * log(accuracy))

cellq2d = CellListDirQ2D(info, r_c, boundary, 1)
interaction_short = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, r_c, n_t)
interaction_long = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, (L, L, 10.0), false, accuracy, α, n_atoms, k_c, 0)

@time Es = QuasiEwald_Es(interaction_short, cellq2d, sys, info)
@time El = QuasiEwald_El(interaction_long, sortz, sys, info)

Es + El