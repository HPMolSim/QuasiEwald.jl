# n_atoms = 1500;

# L = (100.0, 100.0, 30.0);
# L_x, L_y, L_z = L;


# atoms = Vector{Atom{Float64}}()
# for i in 1:1000
#     push!(atoms, Atom(mass = 1.0, charge = 1.0))
# end

# for i in 1:500
#     push!(atoms, Atom(mass = 1.0, charge = -2.0))
# end
# q = [atoms[i].charge for i in 1:n_atoms]

# boundary = Q2dBoundary(L_x, L_y, L_z)
# info = SimulationInfo(n_atoms, atoms, (0.0, L_x, 0.0, L_y, 0.5, L_z - 0.5), boundary; min_r = 2.0, temp = 1.0);
# coords = info.coords;
# z_coords = [coord[3] for coord in coords];
# z_list = sortperm(z_coords);

# γ_1 = 0.95
# γ_2 = 0.95
# ϵ_0 = 1.0


# # rbe off
# accuracy = 1e-4
# α = 1.0
# k_c = sqrt(- 4 * α * log(accuracy))
# r_c = (α * accuracy)^(-1/3)
# n_t = 30
# rbe_p = 250

# intershort = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, false, accuracy, α, n_atoms, r_c, n_t)
# shortfinder = CellListDirQ2D(info, r_c + 1.0, boundary, 100)
# interlong = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, false, accuracy, α, n_atoms, k_c, 0)
# long_finder = SortingFinder(coords)

# mass = [1.0 for i in 1:n_atoms];
# T = Float64
# force_s = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms];
# @time QuasiEwald_Fs!(intershort, shortfinder, atoms, boundary, info.coords, force_s)

# force_l = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms];
# @time QuasiEwald_Fl!(interlong, long_finder, atoms, boundary, info.coords, force_l)


# begin
#     # rbe on
#     accuracy = 1e-3
#     n_t = 30
#     rbe_p = 250

#     α = RBE_α(n_atoms, L_x, L_y, n_t, accuracy, rbe_p)
#     k_c = sqrt(- 4 * α * log(accuracy))
#     r_c = (α * accuracy)^(-1/3)
    

#     intershort = QuasiEwaldShortInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, r_c, n_t)
#     shortfinder = CellListDirQ2D(info, r_c + 1.0, boundary, 100)
#     interlong = QuasiEwaldLongInteraction(γ_1, γ_2, ϵ_0, L, true, accuracy, α, n_atoms, k_c, rbe_p)
#     long_finder = SortingFinder(coords)

#     mass = [1.0 for i in 1:n_atoms];
#     T = Float64
#     force_qem = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
#     @time QuasiEwald_Fs!(intershort, shortfinder, atoms, boundary, info.coords, force_qem)

#     @time QuasiEwald_Fl!(interlong, long_finder, atoms, boundary, info.coords, force_qem)
# end

# # result by ICM
# N_img = 10
# N_real = 10
# sys = IcmSys((γ_2, γ_1), L, N_real, N_img)
# ref_pos, ref_charge = IcmSysInit(sys, coords, q)
# force_jl = IcmForce(sys, coords, q, ref_pos, ref_charge) ./ ϵ_0;



@testset "compare the force with ICM method" begin
    
end