using BenchmarkTools

n_atoms = 1000
q = 2 .* rand(n_atoms) .- 1;
q[end] -= sum(q);

L = (100.0, 100.0, 10.0)

coords = [Point(100.0 * rand(), 100.0 * rand(), 10.0 * rand()) for i in 1:n_atoms]
z_coords = [coord[3] for coord in coords]
z_list = sortperm(z_coords)



γ_1 = 0.8
γ_2 = -0.9
ϵ_0 = 1.0
α = 10.0
accuracy = 1e-4
k_c = sqrt(-4 * α * log(accuracy))

rbe_p = 30
K_set, S = rbe_sampling(L, α, accuracy)

exact_result = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)
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