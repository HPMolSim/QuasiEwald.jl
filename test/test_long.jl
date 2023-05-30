using BenchmarkTools, Plots

n_atoms = 50
q = 2 .* rand(n_atoms) .- 1;
q[end] -= sum(q);

L = (10.0, 10.0, 10.0)

coords = [Point(10.0 * rand(), 10.0 * rand(), 10.0 * rand()) for i in 1:n_atoms]
z_coords = [coord[3] for coord in coords]
z_list = sortperm(z_coords)



γ_1 = 0.3
γ_2 = - 0.4
ϵ_0 = 1.0

accuracy = 1e-4

α = 20.0

α_array = 10.0:5.0:100.0
exact_result_array = []
for α in α_array
    k_c = sqrt(-4 * α * log(accuracy))
    exact_result = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)
    push!(exact_result_array, exact_result - sum(q.^2 .* sqrt(α / π) / 8))
end
plot(α_array, exact_result_array)

energy_ICM = []
N_img_array = 4:2:10
for N_img in N_img_array
    N_real = 20
    ICM_sys = IcmSys((γ_2, γ_1), L, N_real, N_img)
    ref_pos, ref_charge = IcmSysInit(ICM_sys, coords, q)
    energy_jl = IcmEnergy(ICM_sys, coords, q, ref_pos, ref_charge)
    push!(energy_ICM, energy_jl)
end
plot(N_img_array, energy_ICM)

energy_ICM = []
N_real_array = 40:20:200
for N_real in N_real_array
    N_img = 10
    ICM_sys = IcmSys((γ_2, γ_1), L, N_real, N_img)
    ref_pos, ref_charge = IcmSysInit(ICM_sys, coords, q)
    energy_jl = IcmEnergy(ICM_sys, coords, q, ref_pos, ref_charge)
    push!(energy_ICM, energy_jl)
end
plot!(N_real_array, energy_ICM)



rbe_p = 30
K_set, S = rbe_sampling(L, α, accuracy)
sampling_result = []
for i in 1:25
    push!(sampling_result, energy_sum_sampling(q, coords, z_list, L, γ_1, γ_2, ϵ_0, rbe_p, S, K_set))
end
mean(sampling_result)

@benchmark energy_sum_total($q, $coords, $z_list, $L, $γ_1, $γ_2, $ϵ_0, $α, $k_c)
@benchmark energy_sum_sampling($q, $coords, $z_list, $L, $γ_1, $γ_2, $ϵ_0, $rbe_p, $S, $K_set)

k_set = (0.1, 0.1, sqrt(0.02))
container = Container{Float64}(n_atoms)
update_container!(container, k_set, n_atoms, 10.0, coords);
green_element = GreensElement(γ_1, γ_2, 10.0, α);
@btime energy_k_sum($k_set, $q, $coords, $z_list, $green_element, $container)
@btime energy_sum_sampling($q, $coords, $z_list, $L, $γ_1, $γ_2, $ϵ_0, $rbe_p, $S, $K_set)