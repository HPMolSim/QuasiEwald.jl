n_atoms = 100;

L = (100.0, 100.0, 10.0);
L_x, L_y, L_z = L;

ϵ_0 = 1.0

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms
    push!(atoms, Atom(mass = 1.0, charge = (-1.0)^i))
end
q = [(-1.0)^i for i in 1:n_atoms]
boundary = Q2dBoundary(L_x, L_y, L_z)
info = SimulationInfo(n_atoms, atoms, (0.0, L_x, 0.0, L_y, 0.5, L_z - 0.5), boundary; min_r = 2.0, temp = 1.0);
coords = info.coords;
z_coords = [coord[3] for coord in coords];
z_list = sortperm(z_coords);

γ_1 = 0.8
γ_2 = - 0.9
ϵ_0 = 1.0

accuracy = 1e-4
α = 10.0
k_c = sqrt(- 4 * α * log(accuracy))
r_c = 2 * (α * accuracy)^(-1/3)
n_t = 100

intershort = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, false, accuracy, α, n_atoms, r_c, n_t)
shortfinder = CellListDirQ2D(info, r_c + 1.0, boundary, 100)
interlong = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, false, accuracy, α, n_atoms, k_c, 0)
long_finder = SortingFinder(coords)

mass = [1.0 for i in 1:n_atoms];
T = Float64
force_qem = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
@time QuasiEwald_Fs!(intershort, shortfinder, atoms, boundary, info.coords, force_qem)

@time QuasiEwald_Fl!(interlong, long_finder, atoms, boundary, info.coords, force_qem)

# result by ICM
N_img = 10
N_real = 50
sys = IcmSys((γ_2, γ_1), L, N_real, N_img)
ref_pos, ref_charge = IcmSysInit(sys, coords, q)
force_jl = IcmForce(sys, coords, q, ref_pos, ref_charge);

