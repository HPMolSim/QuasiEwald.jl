using Plots, Test, StatsPlots
using ExTinyMD, QuasiEwald

function vecpoint2vec(vecpoint::Vector{Point{3, T}}) where{T}
    vec = Vector{T}()
    for i in 1:length(vecpoint)
        push!(vec, vecpoint[i][1])
        push!(vec, vecpoint[i][2])
        push!(vec, vecpoint[i][3])
    end
    return vec
end

# diff_vec(a, b) = a ⋅ (a - b) / norm(a)
function diff_vec(a::Vector{T}, b::Vector{T}) where{T}
    return dot(a, a .- b) / (norm(a))
end

begin
    n_atoms = 100
    q = 2 .* randn(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms

    Lx, Ly, Lz = (100.0, 100.0, 10.0)
    L = (Lx, Ly, Lz)

    boundary = ExTinyMD.Q2dBoundary(Lx, Ly, Lz)
    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms
        push!(atoms, Atom(mass = 1.0, charge = q[i]))
    end
    info = SimulationInfo(n_atoms, atoms, (0.0, Lx, 0.0, Ly, 0.0, Lz), boundary; min_r = 1.0, temp = 1.0)

    coords = info.coords
    z_coords = [coords[i][3] for i=1:n_atoms]
    z_list = sortperm(z_coords)
    mass = [atoms[i].mass for i in 1:n_atoms]

    accuracy = 1e-4
    α = 0.1
    ϵ_0 = 1.0
    rbe_p = 100
    k_c = sqrt( - 4α * log(accuracy))
    K_set, S = rbe_sampling(L, α, accuracy)

    #benchmark energy non_div
    for (γ_1, γ_2) in [(0.9, 0.9)]
        sort_sampling_array = Vector{Float64}()
        sort_sum = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)

        N_t = 10000
        for i in 1:N_t
            sort_sampling = energy_sum_sampling(q, coords, z_list, L, γ_1, γ_2, ϵ_0, rbe_p, S, K_set)
            push!(sort_sampling_array, sort_sampling)
        end

        round_sort = round(sum(sort_sampling_array) / N_t, digits=2)
        round_exact = round(sort_sum, digits=2)

        sort_sampling_norm = sort_sampling_array ./ sort_sum
        round_sort_norm = round(sum(sort_sampling_norm) / N_t, digits=2)

        var_sort = round(var(sort_sampling_array), digits=2)
        var_sort_norm = round(var(sort_sampling_norm), digits=2)

        plot1 = histogram(sort_sampling_array, bins = 50, label = "sort", xlabel = "value", ylabel = "times", title = "var = $var_sort, mean = $round_sort")
        vline!([sort_sum], linecolor = :red, linestyle = :dash, linewidth = 2, label = "exact = $round_exact")

        plot2 = histogram(sort_sampling_norm, bins = 50, label = "sort", xlabel = "value", ylabel = "times", title = "var_norm = $var_sort_norm")
        vline!([1.0], linecolor = :red, linestyle = :dash, linewidth = 2, label = "1.0")

        plot(plot1, plot2, layout = (1, 2), size = [800, 300], dpi = 500, margin=10Plots.mm)
        savefig("sampling_non_div.pdf")
    end

    #benchmark energy_div
    for (γ_1, γ_2) in [(10.0, 10.0)]
        
        Δk = 2π / Lx
        k_0 = log(γ_1 * γ_2) / (2 * Lz)
        ra = RingAngles(k_0, Lx, Ly, Lz, α, k_c, Δk)

        split_sampling_array = Vector{Float64}()
        sort_sampling_array = Vector{Float64}()

        split_sum = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c, ra)
        sort_sum = energy_sum_total(q, coords, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)

        N_t = 10000
        for i in 1:N_t
            split_sampling = energy_sum_sampling(q, coords, z_list, L, γ_1, γ_2, ϵ_0, rbe_p, S, K_set, ra)
            sort_sampling = energy_sum_sampling(q, coords, z_list, L, γ_1, γ_2, ϵ_0, rbe_p, S, K_set)

            push!(split_sampling_array, split_sampling)
            push!(sort_sampling_array, sort_sampling)
        end

        @test split_sum ≈ sort_sum

        round_split = round(sum(split_sampling_array) / N_t, digits=2)
        round_sort = round(sum(sort_sampling_array) / N_t, digits=2)
        round_exact = round(sort_sum, digits=2)

        var_split = round(var(split_sampling_array), digits=2)
        var_sort = round(var(sort_sampling_array), digits=2)

        plot1 = histogram(split_sampling_array, bins = 50, label = "split", xlabel = "value", ylabel = "times", title = "var = $var_split, mean = $round_split")
        vline!([sort_sum], linecolor = :red, linestyle = :dash, linewidth = 2, label = "exact = $round_exact")

        plot2 = histogram(sort_sampling_array, bins = 50, label = "sort", xlabel = "value", ylabel = "times", title = "var = $var_sort, mean = $round_sort")
        vline!([sort_sum], linecolor = :red, linestyle = :dash, linewidth = 2, label = "exact = $round_exact")

        split_sampling_norm = split_sampling_array ./ sort_sum
        sort_sampling_norm = sort_sampling_array ./ sort_sum

        var_split_norm = round(var(split_sampling_norm), digits=2)
        var_sort_norm = round(var(sort_sampling_norm), digits=2)

        plot3 = histogram(split_sampling_norm, bins = 50, label = "split", xlabel = "value", ylabel = "times", title = "var_norm = $var_split_norm")
        vline!([1.0], linecolor = :red, linestyle = :dash, linewidth = 2, label = "1.0")
        plot4 = histogram(sort_sampling_norm, bins = 50, label = "sort", xlabel = "value", ylabel = "times", title = "var_norm = $var_sort_norm")
        vline!([1.0], linecolor = :red, linestyle = :dash, linewidth = 2, label = "1.0")
        plot(plot1, plot2, plot3, plot4, layout = (2, 2), size = [800, 600], dpi = 1000, margin=10Plots.mm)
        savefig("sampling_div.pdf")
    end

    T = Float64
    #benchmark force non_div
    for (γ_1, γ_2) in [(0.9, 0.9)]
        sort_sampling_array = Vector{Vector{Float64}}()

        sort_sum = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
        force_long_total!(q, mass, coords, sort_sum, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)
        sort_sum_vec = vecpoint2vec(sort_sum)


        N_t = 10000
        for i in 1:N_t
            sort_sampling = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
            force_long_sampling!(q, mass, coords, sort_sampling, z_list, L, γ_1, γ_2, ϵ_0, rbe_p, S, K_set)
            push!(sort_sampling_array, vecpoint2vec(sort_sampling))
        end

        sort_mean = zeros(T, 3 * n_atoms)
        for i in 1:3 * n_atoms
            for j in 1:N_t
                sort_mean[i] += sort_sampling_array[j][i] / N_t
            end
        end

        norm_exact = norm(sort_sum_vec)
        round_exact = round(norm_exact, digits=2)

        error_sort = norm(sort_sum_vec .- sort_mean)
        round_error_sort = round(error_sort, digits=2)

        diff_sort = [diff_vec(sort_sum_vec, sort_sampling_array[i]) for i in 1:N_t]
        diff_sort_mean = diff_vec(sort_sum_vec, sort_mean)
        round_diff_sort_mean = round(diff_sort_mean, digits=2)

        histogram(diff_sort, bins = 50, label = "sort", xlabel = "(ΔF ⋅ F) / |F|", ylabel = "times", margin=10Plots.mm, title = "|F| = $round_exact, |mean(ΔF)| = $round_error_sort", dpi = 1000)
        vline!([diff_sort_mean], linecolor = :red, linestyle = :dash, linewidth = 2, label = "Δ = $round_diff_sort_mean")
        savefig("force_sampling_non_div.pdf")
    end

    #benchmark force div
    for (γ_1, γ_2) in [(10.0, 10.0)]
        k_0 = log(γ_1 * γ_2) / (2 * Lz)
        ra = RingAngles(k_0, Lx, Ly, Lz, α, k_c, π/Lx)

        sort_sampling_array = Vector{Vector{Float64}}()
        split_sampling_array = Vector{Vector{Float64}}()
        sort_sum = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
        force_long_total!(q, mass, coords, sort_sum, z_list, L, γ_1, γ_2, ϵ_0, α, k_c)
        sort_sum_vec = vecpoint2vec(sort_sum)


        N_t = 10000
        for i in 1:N_t
            sort_sampling = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
            force_long_sampling!(q, mass, coords, sort_sampling, z_list, L, γ_1, γ_2, ϵ_0, rbe_p, S, K_set)
            push!(sort_sampling_array, vecpoint2vec(sort_sampling))

            split_sampling = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
            force_long_sampling!(q, mass, coords, split_sampling, z_list, L, γ_1, γ_2, ϵ_0, rbe_p, S, K_set, ra)
            push!(split_sampling_array, vecpoint2vec(split_sampling))
        end

        sort_mean = zeros(T, 3 * n_atoms)
        split_mean = zeros(T, 3 * n_atoms)
        for i in 1:3 * n_atoms
            for j in 1:N_t
                sort_mean[i] += sort_sampling_array[j][i] / N_t
                split_mean[i] += split_sampling_array[j][i] / N_t
            end
        end

        norm_exact = norm(sort_sum_vec)
        round_exact = round(norm_exact, digits=2)

        error_sort = norm(sort_sum_vec .- sort_mean)
        round_error_sort = round(error_sort, digits=2)

        error_split = norm(sort_sum_vec .- split_mean)
        round_error_split = round(error_split, digits=2)

        diff_sort = [diff_vec(sort_sum_vec, sort_sampling_array[i]) for i in 1:N_t]
        diff_sort_mean = diff_vec(sort_sum_vec, sort_mean)
        round_diff_sort_mean = round(diff_sort_mean, digits=2)

        diff_split = [diff_vec(sort_sum_vec, split_sampling_array[i]) for i in 1:N_t]
        diff_split_mean = diff_vec(sort_sum_vec, split_mean)
        round_diff_split_mean = round(diff_split_mean, digits=2)

        plot_1 = histogram(diff_sort, bins = 50, label = "sort", xlabel = "(ΔF ⋅ F) / |F|", ylabel = "times", margin=10Plots.mm, title = "|F| = $round_exact, |mean(ΔF)| = $round_error_sort")
        vline!([diff_sort_mean], linecolor = :red, linestyle = :dash, linewidth = 2, label = "Δ = $round_diff_sort_mean")

        plot_2 = histogram(diff_split, bins = 50, label = "split", xlabel = "(ΔF ⋅ F) / |F|", ylabel = "times", margin=10Plots.mm, title = "|mean(ΔF)| = $round_error_split")
        vline!([diff_split_mean], linecolor = :red, linestyle = :dash, linewidth = 2, label = "Δ = $round_diff_split_mean")

        plot(plot_1, plot_2, layout = (1, 2), size = [800, 300], margin=10Plots.mm, dpi = 1000)
        savefig("force_sampling_div.pdf")
    end
end