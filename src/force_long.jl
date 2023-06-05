export QuasiEwald_Fl!

function QuasiEwald_Fl!(interaction::QuasiEwaldLongInteraction{T, TI}, neighborfinder::SortingFinder{T, TI}, atoms::Vector{ExTinyMD.Atom{T}}, boundary::ExTinyMD.Boundary{T}, coords::Vector{Point{3, T}}, acceleration::Vector{Point{3, T}}) where {T<:Number, TI<:Integer}

    q = [atom.charge for atom in atoms]
    mass = [atom.mass for atom in atoms]

    if interaction.rbe == true
        force_long_sampling!(q, mass, coords, acceleration, neighborfinder.z_list, interaction.L, interaction.γ_1, interaction.γ_2, interaction.ϵ_0, interaction.rbe_p, interaction.S, interaction.K_set)
    else
        force_long_total!(q, mass, coords, acceleration, neighborfinder.z_list, interaction.L, interaction.γ_1, interaction.γ_2, interaction.ϵ_0, interaction.α, interaction.k_c)
    end

    return nothing
end

function force_long_total!(q::Vector{T}, mass::Vector{T}, coords::Vector{Point{3, T}}, acceleration::Vector{Point{3, T}}, z_list::Vector{TI}, L::NTuple{3, T}, γ_1::T, γ_2::T, ϵ_0::T, α::T, k_c::T) where {T<:Number, TI<:Integer}
    n_atoms = size(coords)[1]
    
    L_x, L_y, L_z = L

    acceleration .+= force_k_sum_0(q, z_list) ./ mass ./ (2 * L_x * L_y * ϵ_0)
    
    n_x_max = TI(round(k_c * L_x / T(2) * π, RoundUp))
    n_y_max = TI(round(k_c * L_y / T(2) * π, RoundUp))

    element = GreensElement(γ_1, γ_2, L_z, α)
    container = Container{T}(n_atoms)
    sum_temp = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]

    for n_x in - n_x_max : n_x_max
        for n_y in - n_y_max : n_y_max
            k_x = T(2) * π * n_x / L_x
            k_y = T(2) * π * n_y / L_y
            k = sqrt(k_x^2 + k_y^2)
            k_set = (k_x, k_y, k)
            if k < k_c && k != 0
                update_container!(container, k_set, n_atoms, L_z, coords)
                erase_sum_temp!(sum_temp)
                force_k_sum_1!(k_set, q, z_list, container, sum_temp, element)
                force_k_sum_2!(k_set, q, z_list, container, sum_temp, element)
                force_k_sum_3!(k_set, q, z_list, container, sum_temp, element)
                force_k_sum_4!(k_set, q, z_list, container, sum_temp, element)
                β = γ_1 * γ_2 * exp(- 2 * k * L_z) - 1
                acceleration .+= sum_temp .* (exp(- k*k / (4 * α)) / (2 * L_x * L_y * ϵ_0 * β)) ./ mass
            end
        end
    end

    return nothing
end

function force_long_sampling!(q::Vector{T}, mass::Vector{T}, coords::Vector{Point{3, T}}, acceleration::Vector{Point{3, T}}, z_list::Vector{TI}, L::NTuple{3, T}, γ_1::T, γ_2::T, ϵ_0::T, rbe_p::TI, S::T, K_set::Vector{NTuple{3, T}}) where {T<:Number, TI<:Integer}
    n_atoms = size(coords)[1]
    
    L_x, L_y, L_z = L

    acceleration .+= force_k_sum_0(q, z_list) ./ mass ./ (2 * L_x * L_y * ϵ_0)

    element = GreensElement(γ_1, γ_2, L_z, α)
    container = Container{T}(n_atoms)
    sum_temp = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]

    for i in 1:rbe_p
        k_set = K_set[rand(1:size(K_set)[1])]
        update_container!(container, k_set, n_atoms, L_z, coords)
        erase_sum_temp!(sum_temp)
        force_k_sum_1!(k_set, q, z_list, container, sum_temp, element)
        force_k_sum_2!(k_set, q, z_list, container, sum_temp, element)
        force_k_sum_3!(k_set, q, z_list, container, sum_temp, element)
        force_k_sum_4!(k_set, q, z_list, container, sum_temp, element)
        β = γ_1 * γ_2 * exp(- 2 * k * L_z) - 1
        acceleration .+= sum_temp .* (S / rbe_p / (2 * L_x * L_y * ϵ_0 * β)) ./ mass
    end

    return nothing
end

function force_k_sum_0(q::Vector{T}, z_list::Vector{TI}) where{T <: Number, TI<:Integer}
    n_atoms = length(z_list)

    Q1 = zeros(T, n_atoms)
    Q2 = zeros(T, n_atoms)

    for i in 2:n_atoms
        lf = z_list[i - 1]
        Q1[i] = Q1[i - 1] + q[lf]

        ib = n_atoms - i + 1
        lb = z_list[ib + 1]
        Q2[ib] = Q2[ib + 1] + q[lb]
    end

    sum_k0 = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
    for i in 1:n_atoms
        l = z_list[i]
        sum_k0[l] = Point(zero(T), zero(T), q[l] * (Q1[i] - Q2[i]))
    end
    return sum_k0
end

# to avoid allocation, force_k_sum_1234 will share the same sum_temp structure to avoid allocations

function erase_sum_temp!(sum_temp::Vector{Point{3, T}}) where T
    for i in 1:length(sum_temp)
        sum_temp[i] *= zero(T)
    end
    return nothing
end

function force_k_sum_1!(k_set::NTuple{3, T}, q::Vector{T}, z_list::Vector{TI}, container::Container{T}, sum_temp::Vector{Point{3, T}}, element::GreensElement{T}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)
    k_x, k_y, k = k_set

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

    for i in 1:n_atoms
        l = z_list[i]
        q_l = q[l]
        sum_ri = - q_l * (
            SIN_list[l] * (EXP_list_2[l] * C1[i] + EXP_list_1[l] * C2[i]) -
            COS_list[l] * (EXP_list_2[l] * S1[i] + EXP_list_1[l] * S2[i]) )

        sum_zi = q_l * (
            COS_list[l] * ( - EXP_list_2[l] * C1[i] + EXP_list_1[l] * C2[i]) +
            SIN_list[l] * ( - EXP_list_2[l] * S1[i] + EXP_list_1[l] * S2[i]) )

        sum_temp[l] += Point(k_x * sum_ri / k, k_y * sum_ri / k, sum_zi)
    end
    return nothing
end

function force_k_sum_2!(k_set::NTuple{3, T}, q::Vector{T}, z_list::Vector{TI}, container::Container{T}, sum_temp::Vector{Point{3, T}}, element::GreensElement{T}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)
    k_x, k_y, k = k_set

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

    for i in 1:n_atoms
        l = z_list[i]
        q_l = q[l]
        sum_ri = - q_l * (
            SIN_list[l] * (EXP_list_3[l] * C1[i] + EXP_list_2[l] * C2[i]) -
            COS_list[l] * (EXP_list_3[l] * S1[i] + EXP_list_2[l] * S2[i]))

        sum_zi = q_l * (
            COS_list[l] * (EXP_list_3[l] * C1[i] - EXP_list_2[l] * C2[i]) +
            SIN_list[l] * (EXP_list_3[l] * S1[i] - EXP_list_2[l] * S2[i]))

        sum_temp[l] += element.γ_1 * element.γ_2 * Point(k_x * sum_ri / k, k_y * sum_ri / k, sum_zi)
    end
    return nothing
end

function force_k_sum_3!(k_set::NTuple{3, T}, q::Vector{T}, z_list::Vector{TI}, container::Container{T}, sum_temp::Vector{Point{3, T}}, element::GreensElement{T}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)
    k_x, k_y, k = k_set

    COS_list = container.COS_list
    SIN_list = container.SIN_list
    EXP_list = container.EXP_list_2

    C = zero(T)
    S = zero(T)
    for i in 1:n_atoms
        qi = q[i]
        val = qi * EXP_list[i]
        C += val * COS_list[i]
        S += val * SIN_list[i]
    end

    for i in 1:n_atoms
        val = q[i] * EXP_list[i]
        sum_ri = - val * (SIN_list[i] * C - COS_list[i] * S)
        sum_zi = - val * (COS_list[i] * C + SIN_list[i] * S)

        sum_temp[i] += element.γ_1 * Point(k_x * sum_ri / k, k_y * sum_ri / k, sum_zi)
    end

    return nothing
end

function force_k_sum_4!(k_set::NTuple{3, T}, q::Vector{T}, z_list::Vector{TI}, container::Container{T}, sum_temp::Vector{Point{3, T}}, element::GreensElement{T}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)
    k_x, k_y, k = k_set

    COS_list = container.COS_list
    SIN_list = container.SIN_list
    EXP_list = container.EXP_list_4

    C = zero(T)
    S = zero(T)
    for i in 1:n_atoms
        qi = q[i]
        val = qi * EXP_list[i]
        C += val * COS_list[i]
        S += val * SIN_list[i]
    end

    for i in 1:n_atoms
        val = q[i] * EXP_list[i]
        sum_ri = - val * (SIN_list[i] * C - COS_list[i] * S)
        sum_zi = val * (COS_list[i] * C + SIN_list[i] * S)

        sum_temp[i] += element.γ_2 * Point(k_x * sum_ri / k, k_y * sum_ri / k, sum_zi)
    end

    return nothing
end

function force_direct_sum(k_set::NTuple{3, T}, q::Vector{T}, coords::Vector{Point{3, T}}, L_z::T, γ_1::T, γ_2::T) where T
    
    k_x, k_y, k = k_set
    n_atoms = length(q)

    sum_direct = [[Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms] for j in 1:4]
    sumr_1, sumr_2, sumr_3, sumr_4 = [zeros(T, n_atoms) for i in 1:4]
    sumz_1, sumz_2, sumz_3, sumz_4 = [zeros(T, n_atoms) for i in 1:4]
    for i in 1:n_atoms
        for j in 1:n_atoms
            xi, yi, zi = [coords[i][l] for l in 1:3]
            xj, yj, zj = [coords[j][l] for l in 1:3]
            if j != i
                qs = q[i] * q[j] *  sin(k_x * (xi - xj) + k_y * (yi - yj))
                sumr_1[i] += qs * exp(-k * abs(zi - zj)) 
                sumr_2[i] += γ_1 * qs * exp(-k * (zi + zj)) 
                sumr_3[i] += γ_2 * qs * exp(-k * (2 * L_z - zi - zj))
                sumr_4[i] += γ_1 * γ_2 * qs * exp(-k * (2 * L_z - abs(zi - zj)))
            end
            qc = q[i] * q[j] *  cos(k_x * (xi - xj) + k_y * (yi - yj))
            sumz_1[i] += - sign(zi - zj) * qc * exp(-k * abs(zi - zj)) 
            sumz_2[i] += - γ_1 * qc * exp(-k * (zi + zj)) 
            sumz_3[i] += + γ_2 * qc * exp(-k * (2 * L_z - zi - zj))
            sumz_4[i] += + γ_1 * γ_2 * sign(zi - zj) * qc * exp(-k * (2 * L_z - abs(zi - zj)))
        end
    end
    
    sum_x = - (k_x / k) .* [sumr_1, sumr_2, sumr_3, sumr_4]
    sum_y = - (k_y / k) .* [sumr_1, sumr_2, sumr_3, sumr_4]
    sum_z = [sumz_1, sumz_2, sumz_3, sumz_4]

    for i in 1:4
        for j in 1:n_atoms
            sum_direct[i][j] = Point(sum_x[i][j], sum_y[i][j], sum_z[i][j])
        end
    end

    return sum_direct
end