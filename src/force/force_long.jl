function QuasiEwald_Fl!(interaction::QuasiEwaldLongInteraction{T, TI}, neighborfinder::SortingFinder{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}

    atoms = sys.atoms

    for i in 1:length(interaction.q)
        interaction.q[i] = atoms[info.particle_info[i].id].charge
        interaction.mass[i] = atoms[info.particle_info[i].id].mass
        interaction.coords[i] = info.particle_info[i].position
    end

     erase_vector_of_point!(interaction.acceleration)

    if interaction.rbe == true
        if interaction.k_0 > 0
            force_long_sampling!(interaction.q, interaction.mass, interaction.coords, interaction.acceleration, neighborfinder.z_list, interaction.L, interaction.γ_1, interaction.γ_2, interaction.ϵ_0, interaction.rbe_p, interaction.sum_k, interaction.K_set, interaction.ringangles)
        else
            force_long_sampling!(interaction.q, interaction.mass, interaction.coords, interaction.acceleration, neighborfinder.z_list, interaction.L, interaction.γ_1, interaction.γ_2, interaction.ϵ_0, interaction.rbe_p, interaction.sum_k, interaction.K_set)
        end
    else
        force_long_total!(interaction.q, interaction.mass, interaction.coords, interaction.acceleration, neighborfinder.z_list, interaction.L, interaction.γ_1, interaction.γ_2, interaction.ϵ_0, interaction.α, interaction.k_c)
    end

    for i in 1:length(interaction.acceleration)
        info.particle_info[i].acceleration += interaction.acceleration[i]
    end

    return nothing
end

function force_long_k!(k_set::NTuple{3, T}, q::Vector{T}, z_list::Vector{TI}, container::Container{T}, sum_temp::Vector{Point{3, T}}, element::GreensElement{T}, coords::Vector{Point{3, T}}) where {T <: Number, TI <: Integer}
    force_k_sum_1!(k_set, q, z_list, container, sum_temp, coords)
    force_k_sum_2!(k_set, q, z_list, container, sum_temp, element)
    force_k_sum_3!(k_set, q, z_list, container, sum_temp, element)
    force_k_sum_4!(k_set, q, z_list, container, sum_temp, element)
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
                erase_vector_of_point!(sum_temp)
                force_long_k!(k_set, q, z_list, container, sum_temp, element, coords)
                β = γ_1 * γ_2 * exp(- 2 * k * L_z) - 1
                acceleration .+= sum_temp .* (exp(- k*k / (4 * α)) / (2 * L_x * L_y * ϵ_0 * β)) ./ mass
            end
        end
    end

    return nothing
end

function force_long_total!(q::Vector{T}, mass::Vector{T}, coords::Vector{Point{3, T}}, acceleration::Vector{Point{3, T}}, z_list::Vector{TI}, L::NTuple{3, T}, γ_1::T, γ_2::T, ϵ_0::T, α::T, k_c::T, ringangles::RingAngles{T}) where {T<:Number, TI<:Integer}
    @assert γ_1 * γ_2 ≥ one(T)

    n_atoms = size(coords)[1]
    L_x, L_y, L_z = L
    k_0 = log(γ_1 * γ_2) / (2 * L_z)

    acceleration .+= force_k_sum_0(q, z_list) ./ mass ./ (2 * L_x * L_y * ϵ_0)
    
    n_x_max = TI(round(k_c * L_x / T(2) * π, RoundUp))
    n_y_max = TI(round(k_c * L_y / T(2) * π, RoundUp))

    element = GreensElement(γ_1, γ_2, L_z, α)
    container = Container{T}(n_atoms)
    container_k0 = Container{T}(n_atoms)

    sum_temp = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
    sum_temp_k0 = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]

    for n_x in - n_x_max : n_x_max
        for n_y in - n_y_max : n_y_max
            k_x = T(2) * π * n_x / L_x
            k_y = T(2) * π * n_y / L_y
            k = sqrt(k_x^2 + k_y^2)
            k_set = (k_x, k_y, k)
            kn0_angle = ringangles.ring_angles[nearest_angle_indice(k_x, k_y, ringangles.ring_angles)]
            kn0_set = (k_0 * cos(kn0_angle), k_0 * sin(kn0_angle), k_0)
            if k < k_c && k != 0
                update_container!(container, k_set, n_atoms, L_z, coords)
                update_container!(container_k0, kn0_set, n_atoms, L_z, coords)

                erase_vector_of_point!(sum_temp)
                erase_vector_of_point!(sum_temp_k0)

                force_long_k!(k_set, q, z_list, container, sum_temp, element, coords)
                force_long_k!(kn0_set, q, z_list, container_k0, sum_temp_k0, element, coords)

                sum_temp .-= sum_temp_k0

                β = γ_1 * γ_2 * exp(- 2 * k * L_z) - 1
                acceleration .+= sum_temp .* (exp(- k*k / (4 * α)) / (2 * L_x * L_y * ϵ_0 * β)) ./ mass
            end
        end
    end

    # the divergent part
    for sector_id in 1:length(ringangles.ring_angles) - 1
        kn0_angle = ringangles.ring_angles[sector_id]
        kn0_set = (k_0 * cos(kn0_angle), k_0 * sin(kn0_angle), k_0)
        update_container!(container_k0, kn0_set, n_atoms, L_z, coords)
        erase_vector_of_point!(sum_temp_k0)
        force_long_k!(kn0_set, q, z_list, container_k0, sum_temp_k0, element, coords)

        acceleration .+= sum_temp_k0 .* (ringangles.sectors_sum[sector_id] / (2 * L_x * L_y * ϵ_0)) ./ mass
    end

    return nothing
end

function force_long_sampling!(q::Vector{T}, mass::Vector{T}, coords::Vector{Point{3, T}}, acceleration::Vector{Point{3, T}}, z_list::Vector{TI}, L::NTuple{3, T}, γ_1::T, γ_2::T, ϵ_0::T, rbe_p::TI, S::T, K_set::Vector{NTuple{3, T}}) where {T<:Number, TI<:Integer}
    n_atoms = size(coords)[1]
    
    L_x, L_y, L_z = L

    acceleration .+= force_k_sum_0(q, z_list) ./ mass ./ (2 * L_x * L_y * ϵ_0)

    element = GreensElement(γ_1, γ_2, L_z, one(T))
    container = Container{T}(n_atoms)
    sum_temp = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]

    for i in 1:rbe_p
        k_set = K_set[rand(1:size(K_set)[1])]
        k_x, k_y, k = k_set
        update_container!(container, k_set, n_atoms, L_z, coords)
        erase_vector_of_point!(sum_temp)
        force_long_k!(k_set, q, z_list, container, sum_temp, element, coords)
        β = γ_1 * γ_2 * exp(- 2 * k * L_z) - 1
        acceleration .+= sum_temp .* (S / rbe_p / (2 * L_x * L_y * ϵ_0 * β)) ./ mass
    end

    return nothing
end

function force_long_sampling!(q::Vector{T}, mass::Vector{T}, coords::Vector{Point{3, T}}, acceleration::Vector{Point{3, T}}, z_list::Vector{TI}, L::NTuple{3, T}, γ_1::T, γ_2::T, ϵ_0::T, rbe_p::TI, S::T, K_set::Vector{NTuple{3, T}}, ringangles::RingAngles{T}) where {T<:Number, TI<:Integer}

    n_atoms = size(coords)[1]
    L_x, L_y, L_z = L
    k_0 = log(γ_1 * γ_2) / (2 * L_z)

    acceleration .+= force_k_sum_0(q, z_list) ./ mass ./ (2 * L_x * L_y * ϵ_0)

    element = GreensElement(γ_1, γ_2, L_z, one(T))
    container = Container{T}(n_atoms)
    container_k0 = Container{T}(n_atoms)

    sum_temp = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
    sum_temp_k0 = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]

    for i in 1:rbe_p
        k_set = K_set[rand(1:size(K_set)[1])]
        k_x, k_y, k = k_set
        kn0_angle = ringangles.ring_angles[nearest_angle_indice(k_x, k_y, ringangles.ring_angles)]
        kn0_set = (k_0 * cos(kn0_angle), k_0 * sin(kn0_angle), k_0)

        update_container!(container, k_set, n_atoms, L_z, coords)
        update_container!(container_k0, kn0_set, n_atoms, L_z, coords)

        erase_vector_of_point!(sum_temp)
        erase_vector_of_point!(sum_temp_k0)

        force_long_k!(k_set, q, z_list, container, sum_temp, element, coords)
        force_long_k!(kn0_set, q, z_list, container_k0, sum_temp_k0, element, coords)

        sum_temp .-= sum_temp_k0

        β = γ_1 * γ_2 * exp(- 2 * k * L_z) - 1
        acceleration .+= sum_temp .* (S / rbe_p / (2 * L_x * L_y * ϵ_0 * β)) ./ mass
    end

    # the divergent part
    for sector_id in 1:length(ringangles.ring_angles) - 1
        kn0_angle = ringangles.ring_angles[sector_id]
        kn0_set = (k_0 * cos(kn0_angle), k_0 * sin(kn0_angle), k_0)
        update_container!(container_k0, kn0_set, n_atoms, L_z, coords)
        erase_vector_of_point!(sum_temp_k0)
        force_long_k!(kn0_set, q, z_list, container_k0, sum_temp_k0, element, coords)

        acceleration .+= sum_temp_k0 .* (ringangles.sectors_sum[sector_id] / (2 * L_x * L_y * ϵ_0)) ./ mass
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

function erase_vector_of_point!(vecter_of_point::Vector{Point{3, T}}) where T
    for i in 1:length(vecter_of_point)
        vecter_of_point[i] = Point(zero(T), zero(T), zero(T))
    end
    return nothing
end

function force_k_sum_1!(k_set::NTuple{3, T}, q::Vector{T}, z_list::Vector{TI}, container::Container{T}, sum_temp::Vector{Point{3, T}}, coords::Vector{Point{3, T}}) where {T <: Number, TI <: Integer}
    n_atoms = length(z_list)
    k_x, k_y, k = k_set

    EXP_P_list = container.EXP_P_list
    EXP_N_list = container.EXP_N_list
    A = container.A
    B = container.B

    A[1] = zero(Complex{T})
    j0 = z_list[1]
    A[2] = q[j0] * EXP_N_list[j0]

    for i in 3:n_atoms
        j = z_list[i - 1]
        l = z_list[i - 2]
        zl = coords[l][3]
        zj = coords[j][3]
        A[i] = A[i - 1] * exp(k * (zl - zj)) + q[j] * EXP_N_list[j]
    end

    B[n_atoms] = zero(Complex{T})
    jn = z_list[n_atoms]
    B[n_atoms - 1] = q[jn] * EXP_N_list[jn]

    for i in n_atoms - 2:-1:1
        j = z_list[i + 1]
        l = z_list[i + 2]
        zl = coords[l][3]
        zj = coords[j][3]
        B[i] = B[i + 1] * exp( - k * (zl - zj)) + q[j] * EXP_N_list[j]
    end

    for i in 2:n_atoms
        j = z_list[i]
        l = z_list[i - 1]
        zj = coords[j][3]
        zl = coords[l][3]
        t = q[j] * EXP_P_list[j] * exp( - k * (zj - zl)) * A[i]
        sum_ri = real(1.0im * t)
        sum_zi = real( - t)
        sum_temp[j] += Point(k_x * sum_ri / k, k_y * sum_ri / k, sum_zi)
    end

    for i in 1:n_atoms - 1
        j = z_list[i]
        l = z_list[i + 1]
        zj = coords[j][3]
        zl = coords[l][3]
        t = q[j] * EXP_P_list[j] * exp(k * (zj - zl)) * B[i]
        sum_ri = real(1.0im * t)
        sum_zi = real(t)
        sum_temp[j] += Point(k_x * sum_ri / k, k_y * sum_ri / k, sum_zi) 
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

function force_direct_sum_k0(q::Vector{T}, coords::Vector{Point{3, T}}) where{T <: Number}
    n_atoms = length(q)
    sum_direct = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]

    for i in 1:n_atoms
        for j in 1:n_atoms
            zi = coords[i][3]
            zj = coords[j][3]
            if j != i
                qs = q[i] * q[j]
                sum_direct[i] += Point(zero(T), zero(T), qs * sign(zi - zj))
            end
        end
    end

    return sum_direct
end

function force_direct_sum_k(k_set::NTuple{3, T}, q::Vector{T}, coords::Vector{Point{3, T}}, L_z::T, γ_1::T, γ_2::T) where T
    
    k_x, k_y, k = k_set
    n_atoms = length(q)

    sum_direct = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
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
    
    sum_x = - (k_x / k) .* (sumr_1 + sumr_2 + sumr_3 + sumr_4)
    sum_y = - (k_y / k) .* (sumr_1 + sumr_2 + sumr_3 + sumr_4)
    sum_z = (sumz_1 + sumz_2 + sumz_3 + sumz_4)

    for j in 1:n_atoms
        sum_direct[j] += Point(sum_x[j], sum_y[j], sum_z[j])
    end

    return sum_direct
end

function force_direct_sum_total(q::Vector{T}, mass::Vector{T}, coords::Vector{Point{3, T}}, L::NTuple{3, T}, γ_1::T, γ_2::T, ϵ_0::T, α::T, k_c::T) where {T<:Number}
    n_atoms = size(coords)[1]

    acceleration = [Point(zero(T), zero(T), zero(T)) for i in 1:n_atoms]
    
    L_x, L_y, L_z = L

    acceleration .+= force_direct_sum_k0(q, coords) ./ mass ./ (2 * L_x * L_y * ϵ_0)
    
    n_x_max = Int(round(k_c * L_x / T(2) * π, RoundUp))
    n_y_max = Int(round(k_c * L_y / T(2) * π, RoundUp))

    for n_x in - n_x_max : n_x_max
        for n_y in - n_y_max : n_y_max
            k_x = T(2) * π * n_x / L_x
            k_y = T(2) * π * n_y / L_y
            k = sqrt(k_x^2 + k_y^2)
            k_set = (k_x, k_y, k)
            if k < k_c && k != 0
                sum_k = force_direct_sum_k(k_set, q, coords, L_z, γ_1, γ_2)
                β = γ_1 * γ_2 * exp(- 2 * k * L_z) - 1
                acceleration .+= sum_k .* (exp(- k*k / (4 * α)) / (2 * L_x * L_y * ϵ_0 * β)) ./ mass
            end
        end
    end

    return acceleration
end