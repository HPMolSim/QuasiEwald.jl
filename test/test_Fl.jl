using BenchmarkTools, Test

n_atoms = 50;
q = randn(n_atoms);
q .-= sum(q) / n_atoms;

L = (10.0, 10.0, 10.0);
L_x, L_y, L_z = L;

ϵ_0 = 1.0

coords = [Point(10.0 * rand(), 10.0 * rand(), 10.0 * rand()) for i in 1:n_atoms];
z_coords = [coord[3] for coord in coords];
z_list = sortperm(z_coords);

γ_1 = 0.2
γ_2 = - 0.3
ϵ_0 = 1.0

accuracy = 1e-4
T = Float64
element = GreensElement(γ_1, γ_2, L_z, α);
container = Container{T}(n_atoms);
k_set = (0.1, 0.1, sqrt(0.02));

force_direct = force_direct_sum(k_set, q, coords, L_z, γ_1, γ_2);

update_container!(container, k_set, n_atoms, L_z, coords);

sum_temp = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms];
force_k_sum_1!(k_set, q, z_list, container, sum_temp, element)
@testset begin 
    for i in 1:n_atoms
    @test isapprox(dist2(force_direct[1][i], sum_temp[i]), 0.0; atol = 1e-8)
    end
end

@benchmark erase_sum_temp!($sum_temp)
@benchmark force_k_sum_1!($k_set, $q, $z_list, $container, $sum_temp,$element)

erase_sum_temp!(sum_temp)

force_k_sum_2!(k_set, q, z_list, container, sum_temp, element)
@testset begin 
    for i in 1:n_atoms
    @test isapprox(dist2(force_direct[4][i], sum_temp[i]), 0.0; atol = 1e-8)
    end
end

@benchmark force_k_sum_2!($k_set, $q, $z_list, $container, $sum_temp,$element)

erase_sum_temp!(sum_temp)

force_k_sum_3!(k_set, q, z_list, container, sum_temp, element)
@testset begin 
    for i in 1:n_atoms
    @test isapprox(dist2(force_direct[2][i], sum_temp[i]), 0.0; atol = 1e-8)
    end
end

@benchmark force_k_sum_3!($k_set, $q, $z_list, $container, $sum_temp,$element)

erase_sum_temp!(sum_temp)

force_k_sum_4!(k_set, q, z_list, container, sum_temp, element)
@testset begin 
    for i in 1:n_atoms
    @test isapprox(dist2(force_direct[3][i], sum_temp[i]), 0.0; atol = 1e-8)
    end
end

@benchmark force_k_sum_4!($k_set, $q, $z_list, $container, $sum_temp,$element)

erase_sum_temp!(sum_temp)
force_k_sum_1!(k_set, q, z_list, container, sum_temp, element)
force_k_sum_2!(k_set, q, z_list, container, sum_temp, element)
force_k_sum_3!(k_set, q, z_list, container, sum_temp, element)
force_k_sum_4!(k_set, q, z_list, container, sum_temp, element)

F_d = sum(force_direct)


@testset begin 
    for i in 1:n_atoms
    @test isapprox(dist2(F_d[i], sum_temp[i]), 0.0; atol = 1e-8)
    end
end


mass = [1.0 for i in 1:n_atoms];
force = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
α = 20.0
k_c = 40.0
force_long_total!(q, mass, coords, force, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)

N_img = 10
N_real = 50
sys = IcmSys((γ_2, γ_1), L, N_real, N_img)
ref_pos, ref_charge = IcmSysInit(sys, coords, q)
force_jl = IcmForce(sys, coords, q, ref_pos, ref_charge)