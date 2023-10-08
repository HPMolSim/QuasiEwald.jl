@testset "compare sort sum with direct sum" begin
    n_atoms = 100
    q0 = 2 .* randn(n_atoms) .- 1.0
    q0 = q0 .- sum(q0) / n_atoms
    ϵ_0 = 1.0
    α = 1.0
    k_c = 1.1

    k_set = (0.3, 0.4, 0.5)

    Lx, Ly, Lz = (100.0, 100.0, 10.0)
    L = (Lx, Ly, Lz)

    boundary = ExTinyMD.Q2dBoundary(Lx, Ly, Lz)
    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms
        push!(atoms, Atom(mass = 1.0, charge = q0[i]))
    end
    info = SimulationInfo(n_atoms, atoms, (0.0, Lx, 0.0, Ly, 0.0, Lz), boundary; min_r = 1.0, temp = 1.0)

    coords = [p_info.position for p_info in info.particle_info]
    z_coords = [coords[i][3] for i=1:n_atoms]
    z_list = sortperm(z_coords)

    q = [atoms[p_info.id].charge for p_info in info.particle_info]
    mass = [atoms[p_info.id].mass for p_info in info.particle_info]

    T = Float64

    for (γ_1, γ_2) in [(0.0, 0.0), (0.4, 0.5), (0.4, -0.5), (-0.4, -0.5)]

        element = GreensElement(γ_1, γ_2, Lz, α)
        container = Container{T}(n_atoms)
        sort_sum_force_k = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]

        update_container!(container, k_set, n_atoms, Lz, coords)
        force_long_k!(k_set, q, z_list, container, sort_sum_force_k, element)

        direct_sum_force_k = force_direct_sum_k(k_set, q, coords, Lz, γ_1, γ_2)

        @testset "compart sum_k" begin
            for i in 1:n_atoms
                for j in 1:3
                    @test sort_sum_force_k[i][j] ≈ direct_sum_force_k[i][j]
                end
            end
        end

        @testset "compart sum_k0" begin
            sum_k0 = force_k_sum_0(q, z_list)
            direct_sum_k0 = force_direct_sum_k0(q, coords)
            for i in 1:n_atoms
                for j in 1:3
                    @test sum_k0[i][j] ≈ direct_sum_k0[i][j]
                end
            end
        end
        
        a_sort = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
        force_long_total!(q, mass, coords, a_sort, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)

        a_dir = force_direct_sum_total(q, mass, coords, L, γ_1, γ_2, ϵ_0, α, k_c)

        @testset "compart sum_total" begin
            for i in 1:n_atoms
                for j in 1:3
                    @test a_sort[i][j] ≈ a_dir[i][j]
                end
            end
        end
    end

    for (γ_1, γ_2) in [(10.0, 10.0), (-10.0, -10.0)]
        a_sort = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
        force_long_total!(q, mass, coords, a_sort, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)

        a_dir = force_direct_sum_total(q, mass, coords, L, γ_1, γ_2, ϵ_0, α, k_c)

        k_0 = log(γ_1 * γ_2) / (2 * Lz)
        ra = RingAngles(k_0, Lx, Ly, Lz, α, k_c, π/Lx)
        a_split = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
        force_long_total!(q, mass, coords, a_split, z_list, L, γ_1, γ_2, ϵ_0, α, k_c, ra)

        @testset "compart split" begin
            for i in 1:n_atoms
                for j in 1:3
                    @test a_sort[i][j] ≈ a_dir[i][j]
                    @test a_split[i][j] ≈ a_dir[i][j]
                    @test a_split[i][j] ≈ a_sort[i][j]
                end
            end
        end
    end
end