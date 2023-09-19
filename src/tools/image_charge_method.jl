export IcmSysInit, IcmEnergy, IcmForce

function IcmSysInit(sys::IcmSys, position::Vector{Point{3, T}}, charge::Vector{T}) where T <: Number
    # @assert sum(charge) == 0
    γ_up = sys.γ[1]
    γ_down = sys.γ[2]
    reflect_position = Vector{Point{3, T}}()
    reflect_charge = Vector{T}()

    for i = 1:length(charge)
        x, y, z = position[i]
        q = charge[i]
        position_up = Vector{T}()
        position_down = Vector{T}()
        charge_up = Vector{T}()
        charge_down = Vector{T}()
        push!(position_up, 2.0 * sys.L[3] - z)
        push!(position_down, - z)
        push!(charge_up, γ_up * q)
        push!(charge_down, γ_down * q)
        for m in 2:sys.N_img
            push!(position_up, 2 * sys.L[3] - position_down[m-1])
            push!(position_down, - position_up[m-1])
            push!(charge_up, γ_up * charge_down[m-1])
            push!(charge_down, γ_down * charge_up[m-1])
        end
        push!(reflect_position, Point(x, y, z))
        push!(reflect_charge, q)
        for m in 1:sys.N_img
            push!(reflect_position, Point(x, y, position_up[m]))
            push!(reflect_charge, charge_up[m])
            push!(reflect_position, Point(x, y, position_down[m]))
            push!(reflect_charge, charge_down[m])
        end
    end

    return reflect_position, reflect_charge
end

function IcmEnergy(sys::IcmSys, position::Vector{Point{3, T}}, charge::Vector{T}, reflect_position::Vector{Point{3, T}}, reflect_charge::Vector{T}) where T<:Number
    energy = zero(T)
    for i in 1:length(charge)
        q_i = charge[i]
        pos_i = position[i]
        for j in 1:length(reflect_charge)
            q_j = reflect_charge[j]
            pos_j = reflect_position[j]
            for mx in -sys.N_real:sys.N_real
                for my in -sys.N_real:sys.N_real
                    r = dist2(pos_i, pos_j + Point(mx * sys.L[1], my * sys.L[2], 0.0))
                    if r != 0.0
                        energy += q_i * q_j / (8π * sqrt(r))
                    end
                end
            end
        end
    end
    return energy
end

function IcmForce(sys::IcmSys, position::Vector{Point{3, T}}, charge::Vector{T}, reflect_position::Vector{Point{3, T}}, reflect_charge::Vector{T}) where T<:Number
    force = [Point(zero(T), zero(T), zero(T)) for _=1:length(charge)]
    for i in 1:length(charge)
        q_i = charge[i]
        pos_i = position[i]
        for j in 1:length(reflect_charge)
            q_j = reflect_charge[j]
            pos_j = reflect_position[j]
            for mx in -sys.N_real:sys.N_real
                for my in -sys.N_real:sys.N_real
                    force[i] += CoulumbForce(q_i, q_j, pos_i, pos_j + Point(mx * sys.L[1], my * sys.L[2], 0.0))
                end
            end
        end
    end
    return force 
end

function IcmForce_self(sys::IcmSys, position::Vector{Point{3, T}}, charge::Vector{T}, reflect_position::Vector{Point{3, T}}, reflect_charge::Vector{T}) where T<:Number
    force_self = [Point(zero(T), zero(T), zero(T)) for _=1:length(charge)]
    for i in 1:length(charge)
        q_i = charge[i]
        pos_i = position[i]
        for j in 1:length(reflect_charge)
            q_j = reflect_charge[j]
            pos_j = reflect_position[j]
            if pos_j[1] ≈ pos_i[1] && pos_j[2] ≈ pos_i[2]
                for mx in -sys.N_real:sys.N_real
                    for my in -sys.N_real:sys.N_real
                        force_self[i] += CoulumbForce(q_i, q_j, pos_i, pos_j + Point(mx * sys.L[1], my * sys.L[2], 0.0))
                    end
                end
            end
        end
    end
    return force_self 
end

function CoulumbForce(q_i::T, q_j::T, coo_i::Point{3, T}, coo_j::Point{3, T}) where T<:Number
    r = dist2(coo_i, coo_j)
    if iszero(r) == false
        rho = sqrt(r)
        F = q_i * q_j / (4π * r)
        angle = (coo_i - coo_j) * (1/rho)
        return F * angle
    else
        return Point(0.0, 0.0, 0.0)
    end
end