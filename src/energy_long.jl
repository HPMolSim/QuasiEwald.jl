function energy_k_sum(k_set::NTuple{3, T}, q::Vector{T}, coords::Vector{Point{3, T}}, z_list::Vector{TI}, element::GreensElement{T}) where{T <: Number, TI<:Integer}
    k_x, k_y, k = k_set
    n_atoms = size(coords)[1]
    L_z = element.L_z
    # q = [atoms.charge[i] for i in 1:n_atoms]
    α = element.α

    γ_1 = element.γ_1
    γ_2 = element.γ_2
    β = 1 / (γ_1 * γ_2 * exp(- 2 * k * L_z) - 1)

    # there are two terms to be summed
    # sum_1 = beta(k) exp(-k z_n)
    # sum_2 = beta(k) g_1 g_2 exp(-k(2 * L_z - z_n))

    # first for sum_1
    # C1 for the forward process and C2 for the backward, same for S1/2
    # C1[i] = \sum_{j<i} q_j exp(+k z_j) cos(k \rho_j)
    # C2[i] = \sum_{j>i} q_j exp(-k z_j) cos(k \rho_j)
    # S1[i] = \sum_{j<i} q_j exp(+k z_j) sin(k \rho_j)
    # S2[i] = \sum_{j>i} q_j exp(-k z_j) sin(k \rho_j)
    C1 = zeros(T, n_atoms)
    C2 = zeros(T, n_atoms)
    S1 = zeros(T, n_atoms)
    S2 = zeros(T, n_atoms)

    for i in 1:n_atoms-1
        #forward process
        lf = z_list[i]
        q_lf = q[lf]
        x_lf, y_lf, z_lf = coords[lf]
        forward_val = q_lf * exp(k * z_lf)
        C1[i + 1] = C1[i] + forward_val * cos(k_x * x_lf + k_y * y_lf)
        S1[i + 1] = S1[i] + forward_val * sin(k_x * x_lf + k_y * y_lf)
    
        #backward process
        back_i = n_atoms - i
        lb = z_list[back_i + 1]
        q_lb = q[lb]
        x_lb, y_lb, z_lb = coords[lb]
        backward_val = q_lb * exp(-k * z_lb)
        C2[back_i] = C2[back_i + 1] + backward_val * cos(k_x * x_lb + k_y * y_lb)
        S2[back_i] = S2[back_i + 1] + backward_val * sin(k_x * x_lb + k_y * y_lb)
    end

    sum_1 = zero(T)
    # sum = \sum_i q_i cos(k \rho_i)(exp(k z_i) C1[i] + exp(-k z_i) C2[i]) +
    #              q_i sin(k \rho_i)(exp(k z_i) S1[i] + exp(-k z_i) S2[i])
    for i in 1:n_atoms
        l = z_list[i]
        q_l = q[l]
        x_l, y_l, z_l = coords[l]
        sum_1 += q_l * (
            cos(k_x * x_l + k_y * y_l) * (
                exp(-k * z_l) * C1[i] + exp(+k * z_l) * C2[i]
            ) + 
            sin(k_x * x_l + k_y * y_l) * (
                exp(-k * z_l) * S1[i] + exp(+k * z_l) * S2[i]
            ) + 
            q_l
        )
    end

    # here we will handle sum_2
    # sum_2 = sum_i q_i (sum_{j<i} q_j exp(-k(2*L_z - (z_i - z_j))) + sum_{j>i} q_j exp(-k(2*L_z - (z_j - z_i))))
    # = sum_i q_i (exp(-k(2*L_z - z_i)) sum_{j<i} q_j exp(-k * z_j) + exp(-k * z_i) sum_{j>i} q_j exp(-k(2*L_z - z_j)))
    # C1 for the forward process and C2 for the backward, same for S1/2
    # C1[i] = \sum_{j<i} q_j exp(-k z_j) cos(k \rho_j)
    # C2[i] = \sum_{j>i} q_j exp(-k (2L_z - z_j)) cos(k \rho_j)
    # S1[i] = \sum_{j<i} q_j exp(-k z_j) sin(k \rho_j)
    # S2[i] = \sum_{j>i} q_j exp(-k (2L_z - z_j)) sin(k \rho_j)
    C1[1] = zero(T)
    S1[1] = zero(T)
    C2[end] = zero(T)
    S2[end] = zero(T)

    for i in 1:n_atoms-1
        #forward process
        lf = z_list[i]
        q_lf = q[lf]
        x_lf, y_lf, z_lf = coords[lf]
        forward_val = q_lf * exp(- k * z_lf)
        C1[i + 1] = C1[i] + forward_val * cos(k_x * x_lf + k_y * y_lf)
        S1[i + 1] = S1[i] + forward_val * sin(k_x * x_lf + k_y * y_lf)
    
        #backward process
        back_i = n_atoms - i
        lb = z_list[back_i + 1]
        q_lb = q[lb]
        x_lb, y_lb, z_lb = coords[lb]
        backward_val = q_lb * exp(-k * (2 * L_z - z_lb))
        C2[back_i] = C2[back_i + 1] + backward_val * cos(k_x * x_lb + k_y * y_lb)
        S2[back_i] = S2[back_i + 1] + backward_val * sin(k_x * x_lb + k_y * y_lb)
    end

    sum_2 = zero(T)
    # sum = \sum_i q_i cos(k \rho_i)(exp(k z_i) C1[i] + exp(-k z_i) C2[i]) +
    #              q_i sin(k \rho_i)(exp(k z_i) S1[i] + exp(-k z_i) S2[i])
    for i in 1:n_atoms
        l = z_list[i]
        q_l = q[l]
        x_l, y_l, z_l = coords[l]
        sum_2 += q_l * (
            cos(k_x * x_l + k_y * y_l) * (
                exp(-k * (2 * L_z - z_l)) * C1[i] + exp(-k * z_l) * C2[i]
            ) + 
            sin(k_x * x_l + k_y * y_l) * (
                exp(-k * (2 * L_z - z_l)) * S1[i] + exp(-k * z_l) * S2[i]
            ) + 
            q_l * exp(-2 * k * L_z)
        )
    end

    C = zero(T)
    S = zero(T)
    for i in 1:n_atoms
        qi = q[i]
        x_i, y_i, z_i = coords[i]
        val = qi * exp(- k * z_i)
        C += val * cos(k_x * x_i + k_y * y_i)
        S += val * sin(k_x * x_i + k_y * y_i)
    end

    sum_3 = zero(T)
    for i in 1:n_atoms
        qi = q[i]
        x_i, y_i, z_i = coords[i]
        val = qi * exp(- k * z_i)
        ci = cos(k_x * x_i + k_y * y_i)
        si = sin(k_x * x_i + k_y * y_i)
        sum_3 += val * (ci * C + si * S)
    end

    C = zero(T)
    S = zero(T)
    for i in 1:n_atoms
        qi = q[i]
        x_i, y_i, z_i = coords[i]
        val = qi * exp(- k * (L_z - z_i))
        C += val * cos(k_x * x_i + k_y * y_i)
        S += val * sin(k_x * x_i + k_y * y_i)
    end

    sum_4 = zero(T)
    for i in 1:n_atoms
        qi = q[i]
        x_i, y_i, z_i = coords[i]
        val = qi * exp(- k * (L_z - z_i))
        ci = cos(k_x * x_i + k_y * y_i)
        si = sin(k_x * x_i + k_y * y_i)
        sum_4 += val * (ci * C + si * S)
    end

    return (sum_1 + γ_1 * γ_2 * sum_2 + γ_1 * sum_3 + γ_2 * sum_4) * β / k
end


function energy_sum_total(q, coords, α, L_x, L_y, L_z, g_1, g_2, eps_0, k_c)
    n_atoms = size(coords)[1]
    z_coords = [coords[i][3] for i in 1:n_atoms]
    z_list = sortperm(z_coords)

    sum_k0 = energy_k0_sum(q, coords, z_list)
    
    sum_total = 0
    n_x_max = trunc(Int, k_c * L_x / 2 * π)
    n_y_max = trunc(Int, k_c * L_y / 2 * π)

    green_element = greens_element_init(g_1, g_2, L_z, α)

    for n_x in - n_x_max : n_x_max
        for n_y in - n_y_max : n_y_max
            k_x = 2 * π * n_x / L_x
            k_y = 2 * π * n_y / L_y
            k = sqrt(k_x^2 + k_y^2)
            k_set = (k_x, k_y, k)
            if k < k_c && k != 0
                sum_k = energy_k_sum(k_set, q, coords, z_list, green_element)
                sum_total += sum_k * exp(-k^2 / (4 * α))
            end
        end
    end

    return - sum_k0 / (4 * L_x * L_y * eps_0) - sum_total / (4 * L_x * L_y * eps_0)
end