export energy_sum_total, energy_sum_sampling, QuasiEwald_El

function QuasiEwald_El(interaction::QuasiEwaldLongInteraction{T, TI}, neighbor::SortingFinder{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    update_finder!(neighbor, info)
    
    q = [atom.charge for atom in sys.atoms]

    if interaction.rbe == true
        return energy_sum_sampling(q, info.coords, neighbor.z_list, interaction.L, interaction.γ_1, interaction.γ_2, interaction.ϵ_0, interaction.rbe_p, interaction.S, interaction.K_set)
    else
        return energy_sum_total(q, info.coords, neighbor.z_list, interaction.L, interaction.γ_1, interaction.γ_2, interaction.ϵ_0, interaction.α, interaction.k_c)
    end
end

@inbounds function energy_k_sum_0(q::Vector{T}, coords::Vector{Point{3, T}}, z_list::Vector{TI}) where{T <: Number, TI<:Integer}
    n_atoms = length(z_list)

    Q_1 = zeros(T, n_atoms)
    Q_2 = zeros(T, n_atoms)

    for i in 2:n_atoms
        l = z_list[i - 1]
        Q_1[i] = Q_1[i - 1] + q[l]
        Q_2[i] = Q_2[i - 1] + q[l] * coords[l][3]
    end

    sum_k_0 = zero(T)
    for i in 1:n_atoms
        l = z_list[i]
        q_i = q[l]
        z_i = coords[l][3]
        sum_k_0 += q_i * z_i * Q_1[i] - q_i * Q_2[i]
    end

    return 2 * sum_k_0
end

@inbounds function energy_k_sum(k_set::NTuple{3, T}, q::Vector{T}, coords::Vector{Point{3, T}}, z_list::Vector{TI}, element::GreensElement{T}, container::Container{T}) where{T <: Number, TI<:Integer}
    k_x, k_y, k = k_set
    n_atoms = size(coords)[1]
    L_z = element.L_z
    # q = [atoms.charge[i] for i in 1:n_atoms]
    α = element.α

    γ_1 = element.γ_1
    γ_2 = element.γ_2
    β = 1 / (γ_1 * γ_2 * exp(- 2 * k * L_z) - 1)


    sum_1 = energy_k_sum_1(q, z_list, container)
    sum_2 = energy_k_sum_2(q, z_list,  exp(-2*k*L_z), container)
    sum_3 = energy_k_sum_3(q, z_list, container)
    sum_4 = energy_k_sum_4(q, z_list, container)

    return (sum_1 + γ_1 * γ_2 * sum_2 + γ_1 * sum_3 + γ_2 * sum_4) * β / k
end

@inbounds function energy_k_sum_1(q::Vector{T}, z_list::Vector{TI}, container::Container{T}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)

    COS_list = container.COS_list
    SIN_list = container.SIN_list
    EXP_list_1 = container.EXP_list_1
    EXP_list_2 = container.EXP_list_2

    C1 = container.C1
    C2 = container.C2
    S1 = container.S1
    S2 = container.S2

    C1[1] = zero(T)
    S1[1] = zero(T)
    C2[end] = zero(T)
    S2[end] = zero(T)

    #sum_1
    sum_1 = zero(T)
    for i in 1:n_atoms-1
        #forward process
        lf = z_list[i]
        q_lf = q[lf]
        forward_val = q_lf * EXP_list_1[lf]
        C1[i + 1] = C1[i] + forward_val * COS_list[lf]
        S1[i + 1] = S1[i] + forward_val * SIN_list[lf]
    
        #backward process
        back_i = n_atoms - i
        lb = z_list[back_i + 1]
        q_lb = q[lb]
        backward_val = q_lb * EXP_list_2[lb]
        C2[back_i] = C2[back_i + 1] + backward_val * COS_list[lb]
        S2[back_i] = S2[back_i + 1] + backward_val * SIN_list[lb]
    end

    # sum = \sum_i q_i cos(k \rho_i)(exp(k z_i) C1[i] + exp(-k z_i) C2[i]) +
    #              q_i sin(k \rho_i)(exp(k z_i) S1[i] + exp(-k z_i) S2[i])
    for i in 1:n_atoms
        l = z_list[i]
        q_l = q[l]
        sum_1 += q_l * (
            COS_list[l] * (EXP_list_2[l] * C1[i] + EXP_list_1[l] * C2[i]) + 
            SIN_list[l] * (EXP_list_2[l] * S1[i] + EXP_list_1[l] * S2[i]) + 
            q_l
        )
    end

    return sum_1
end

@inbounds function energy_k_sum_2(q::Vector{T}, z_list::Vector{TI}, EXP_Lz::T, container::Container{T}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)

    COS_list = container.COS_list
    SIN_list = container.SIN_list
    EXP_list_3 = container.EXP_list_3
    EXP_list_2 = container.EXP_list_2

    C1 = container.C1
    C2 = container.C2
    S1 = container.S1
    S2 = container.S2

    C1[1] = zero(T)
    S1[1] = zero(T)
    C2[end] = zero(T)
    S2[end] = zero(T)
    
    sum_2 = zero(T)
    
    #sum_2
    for i in 1:n_atoms-1
        #forward process
        lf = z_list[i]
        q_lf = q[lf]
        forward_val = q_lf * EXP_list_2[lf]
        C1[i + 1] = C1[i] + forward_val * COS_list[lf]
        S1[i + 1] = S1[i] + forward_val * SIN_list[lf]
    
        #backward process
        back_i = n_atoms - i
        lb = z_list[back_i + 1]
        q_lb = q[lb]
        backward_val = q_lb * EXP_list_3[lb]
        C2[back_i] = C2[back_i + 1] + backward_val * COS_list[lb]
        S2[back_i] = S2[back_i + 1] + backward_val * SIN_list[lb]
    end

    # sum = \sum_i q_i cos(k \rho_i)(exp(k z_i) C1[i] + exp(-k z_i) C2[i]) +
    #              q_i sin(k \rho_i)(exp(k z_i) S1[i] + exp(-k z_i) S2[i])
    for i in 1:n_atoms
        l = z_list[i]
        q_l = q[l]
        sum_2 += q_l * (
            COS_list[l] * (
                EXP_list_3[l] * C1[i] + EXP_list_2[l] * C2[i]
            ) + 
            SIN_list[l] * (
                EXP_list_3[l] * S1[i] + EXP_list_2[l] * S2[i]
            ) + 
            q_l * EXP_Lz
        )
    end

    return sum_2
end

@inbounds function energy_k_sum_3(q::Vector{T}, z_list::Vector{TI}, container::Container{T}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)

    COS_list = container.COS_list
    SIN_list = container.SIN_list
    EXP_list = container.EXP_list_2

    C = zero(T)
    S = zero(T)
    sum_3 = zero(T)

    for i in 1:n_atoms
        qi = q[i]
        val = qi * EXP_list[i]
        C += val * COS_list[i]
        S += val * SIN_list[i]
    end

    for i in 1:n_atoms
        qi = q[i]
        val = qi * EXP_list[i]
        ci = COS_list[i]
        si = SIN_list[i]
        sum_3 += val * (ci * C + si * S)
    end

    return sum_3
end

@inbounds function energy_k_sum_4(q::Vector{T}, z_list::Vector{TI}, container::Container{T}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)

    COS_list = container.COS_list
    SIN_list = container.SIN_list
    EXP_list = container.EXP_list_4

    C = zero(T)
    S = zero(T)
    sum_4 = zero(T)

    for i in 1:n_atoms
        qi = q[i]
        val = qi * EXP_list[i]
        C += val * COS_list[i]
        S += val * SIN_list[i]
    end

    for i in 1:n_atoms
        qi = q[i]
        val = qi * EXP_list[i]
        ci = COS_list[i]
        si = SIN_list[i]
        sum_4 += val * (ci * C + si * S)
    end

    return sum_4
end


function energy_sum_total(q::Vector{T}, coords::Vector{Point{3, T}}, z_list::Vector{TI}, L::NTuple{3, T}, γ_1::T, γ_2::T, ϵ_0::T, α::T, k_c::T) where {T<:Number, TI<:Integer}
    n_atoms = size(coords)[1]
    L_x, L_y, L_z = L

    sum_k0 = energy_k_sum_0(q, coords, z_list)
    
    sum_total = zero(T)
    n_x_max = TI(round(k_c * L_x / T(2) * π, RoundUp))
    n_y_max = TI(round(k_c * L_y / T(2) * π, RoundUp))

    green_element = GreensElement(γ_1, γ_2, L_z, α)
    container = Container{T}(n_atoms)

    for n_x in - n_x_max : n_x_max
        for n_y in - n_y_max : n_y_max
            k_x = T(2) * π * n_x / L_x
            k_y = T(2) * π * n_y / L_y
            k = sqrt(k_x^2 + k_y^2)
            k_set = (k_x, k_y, k)
            if k < k_c && k != 0
                update_container!(container, k_set, n_atoms, L_z, coords)
                sum_k = energy_k_sum(k_set, q, coords, z_list, green_element, container)
                sum_total += sum_k * exp(- k*k / (4 * α))
            end
        end
    end

    return - (sum_k0 + sum_total) / (T(4) * L_x * L_y * ϵ_0)
end

function energy_sum_sampling(q::Vector{T}, coords::Vector{Point{3, T}}, z_list::Vector{TI}, L::NTuple{3, T}, γ_1::T, γ_2::T, ϵ_0::T, rbe_p::TI, S::T, K_set::Vector{NTuple{3, T}}) where {T<:Number, TI<:Integer}
    n_atoms = size(coords)[1]
    L_x, L_y, L_z = L

    green_element = GreensElement(γ_1, γ_2, L_z, α)
    container = Container{T}(n_atoms)

    sum_k0 = energy_k_sum_0(q, coords, z_list)
    sum_total = zero(T)

    for i in 1:rbe_p
        k_set = K_set[rand(1:size(K_set)[1])]
        sum_k = energy_k_sum(k_set, q, coords, z_list, green_element, container)
        sum_total += sum_k * S / rbe_p
    end

    return - (sum_k0 + sum_total) / (T(4) * L_x * L_y * ϵ_0)
end


# this are three function used to verify our summation method
# they do the summation directly instead of by sorting
function direct_sum_total(q::Vector{T}, coords::Vector{Point{3, T}}, L_x::T, L_y::T, L_z::T, k_c::T, α::T, γ_1::T, γ_2::T) where {T}
    n_atoms = size(coords)[1]

    sum_k0 = zero(T)
    for i in 1:n_atoms
        for j in 1:n_atoms
            sum_k0 += q[i] * q[j] * abs(coords[i][3] - coords[j][3])
        end
    end
    
    sum_k = zero(T)
    n_x_max = trunc(Int, k_c * L_x / 2 * π)
    n_y_max = trunc(Int, k_c * L_y / 2 * π)

    for n_x in - n_x_max : n_x_max
        for n_y in - n_y_max : n_y_max
            k_x = 2 * π * n_x / L_x
            k_y = 2 * π * n_y / L_y
            k = sqrt(k_x^2 + k_y^2)
            if k < k_c && k != 0
                beta = (g_1 * g_2 * exp(-2 * k * L_z) - 1)
                sum_1 = 0
                sum_2 = 0
                sum_3 = 0
                sum_4 = 0
                for i in 1:n_atoms
                    for j in 1:n_atoms
                        xi, yi, zi = coords[i]
                        xj, yj, zj = coords[j]
                        qc = q[i] * q[j] *  cos(k_x * (xi - xj) + k_y * (yi - yj)) / (beta * k)
                        sum_1 += qc * exp(-k * abs(zi - zj)) 
                        sum_2 += g_1 * qc * exp(-k * (zi + zj)) 
                        sum_3 += g_2 * qc * exp(-k * (2 * L_z - zi - zj))
                        sum_4 += g_1 * g_2 * qc * exp(-k * (2 * L_z - abs(zi - zj)))
                    end
                end
                sum_k += (sum_1 + sum_2 + sum_3 + sum_4) * exp(-k^2 / (4 * alpha))
            end
        end
    end

    return -(sum_k0 + sum_k) / (4 * L_x * L_y)
end

function direct_sum_k(k_set::NTuple{3, T}, q::Vector{T}, coords::Vector{Point{3, T}}, element::GreensElement{T}) where {T<:Number}
    n_atoms = size(coords)[1]
    γ_1 = element.γ_1
    γ_2 = element.γ_2
    L_z = element.L_z
    
    k_x, k_y, k = k_set

    β = (γ_1 * γ_2 * exp(-2 * k * L_z) - 1)
    sum_1 = zero(T)
    sum_2 = zero(T)
    sum_3 = zero(T)
    sum_4 = zero(T)
    for i in 1:n_atoms
        for j in 1:n_atoms
            xi, yi, zi = coords[i]
            xj, yj, zj = coords[j]
            qc = q[i] * q[j] *  cos(k_x * (xi - xj) + k_y * (yi - yj)) / (β * k)
            sum_1 += qc * exp(-k * abs(zi - zj)) 
            sum_2 += γ_1 * qc * exp(-k * (zi + zj)) 
            sum_3 += γ_2 * qc * exp(-k * (2 * L_z - zi - zj))
            sum_4 += γ_1 * γ_2 * qc * exp(-k * (2 * L_z - abs(zi - zj)))
        end
    end
    sum_k = (sum_1 + sum_2 + sum_3 + sum_4)

    return sum_k
end

function direct_sum_k_0(q::Vector{T}, coords::Vector{Point{3, T}}) where{T}
    n_atoms = size(coords)[1]

    # k = 0 part
    sum_k0 = zero(T)
    for i in 1:n_atoms
        for j in 1:n_atoms
            sum_k0 += q[i] * q[j] * abs(coords[i][3] - coords[j][3])
        end
    end

    return sum_k0
end