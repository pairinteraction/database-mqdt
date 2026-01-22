import LinearAlgebra
import MQDT
include("angular.jl")


const START_ID = 0

macro timelog(expr)
    quote
        local result, elapsed_time, bytes_allocated, gc_time, memory_counters =
            @timed $(esc(expr))
        @info "$(round(elapsed_time, digits=6)) seconds (allocations: $(Base.format_bytes(bytes_allocated)), $(round(100 * gc_time / elapsed_time, digits=2))% gc time)"
        result
    end
end


function get_relevant_lr(state::MQDT.BasisState)
    inds = findall(state.model.core)
    return state.lr_list[inds]
end

function get_relevant_nu(state::MQDT.BasisState)
    inds = findall(state.model.core)
    return state.nu_list[inds]
end

function all_matrix_element(B::MQDT.BasisArray, parameters::MQDT.Parameters)
    """Calculate all relevant matrix elements for a given basis array B.

    This means dipole, quadrupole, magnetic, and diamagnetic matrix elements.
    """
    k_angular_max = 2  # 2 for now, since we dont calculate octupole matrix elements

    row_col_value = Dict(
        "matrix_elements_d" => Tuple{Int64,Int64,Float64}[],
        "matrix_elements_q" => Tuple{Int64,Int64,Float64}[],
        "matrix_elements_mu" => Tuple{Int64,Int64,Float64}[],
        "matrix_elements_q0" => Tuple{Int64,Int64,Float64}[],
    )

    states_indexed = [(ids, state) for (ids, state) in enumerate(B.states)]
    states_lr = [get_relevant_lr(s) for (_, s) in states_indexed]
    states_nu = [get_relevant_nu(s) for (_, s) in states_indexed]

    states_sorted = sort(
        states_indexed,
        by = x ->
            (minimum(get_relevant_lr(x[2])), minimum(get_relevant_nu(x[2])), x[1]),
    )

    for (i1, (id1, s1)) in enumerate(states_sorted)
        lr1 = states_lr[id1]
        nus1 = states_nu[id1]
        for (id2, s2) in states_sorted[i1:end]
            lr2 = states_lr[id2]
            # Skip if all contributions of the two states are far apart in angular momentum
            if minimum(lr2) - maximum(lr1) > k_angular_max
                continue
            end

            nus2 = states_nu[id2]
            # Skip if all contributions of the two states are far apart in n and None of them is low-n
            if all(abs(nu1-nu2) >= 11 for nu1 in nus1 for nu2 in nus2) &&
               all(nu1 > 25 for nu1 in nus1) &&
               all(nu2 > 25 for nu2 in nus2)
                continue
            end

            m = MQDT.multipole_moments(s1, s2, parameters)
            # multipole_moments returns the matrix elements in the following order
            # electric dipole, electric quadrupole, diamagnetic, magnetic
            table_keys = [
                "matrix_elements_d",
                "matrix_elements_q",
                "matrix_elements_q0",
                "matrix_elements_mu",
            ]
            prefactor_transposed = (-1)^(s2.f - s1.f)

            for (i, key) in enumerate(table_keys)
                if m[i] != 0
                    # start IDs from 0 for consistency with python
                    _id1 = id1 - 1 + START_ID
                    _id2 = id2 - 1 + START_ID
                    push!(row_col_value[key], (_id1, _id2, m[i]))
                    if id1 != id2
                        push!(row_col_value[key], (_id2, _id1, m[i] * prefactor_transposed))
                    end
                end
            end

        end
    end

    return row_col_value
end


function rcv_to_df(row_col_value::Vector{Tuple{Int64,Int64,Float64}})
    """Convert a row_col_value to a DataFrame."""
    id_initial = [m[1] for m in row_col_value]
    id_final = [m[2] for m in row_col_value]
    val = [m[3] for m in row_col_value]
    df = DataFrame(id_initial = id_initial, id_final = id_final, val = val)
    return df
end


function basis_to_df(T::MQDT.BasisArray, P::MQDT.Parameters)
    df = DataFrame(
        id = collect(START_ID:(size(T)-1+START_ID)),
        energy = MQDT.get_e(T, P) / 219474.6313632,  # convert 1/cm to atomic units
        parity = MQDT.get_p(T),
        n = get_n(T, P),
        nu = MQDT.get_nu(T),
        f = MQDT.get_f(T),
        exp_nui = exp_nui(T),
        exp_l = calc_exp_qn(T, "l_tot"),
        exp_j = calc_exp_qn(T, "j_tot"),
        exp_s = calc_exp_qn(T, "s_tot"),
        exp_l_ryd = calc_exp_qn(T, "l_r"),
        exp_j_ryd = calc_exp_qn(T, "j_r"),
        std_nui = std_nui(T),
        std_l = calc_std_qn(T, "l_tot"),
        std_j = calc_std_qn(T, "j_tot"),
        std_s = calc_std_qn(T, "s_tot"),
        std_l_ryd = calc_std_qn(T, "l_r"),
        std_j_ryd = calc_std_qn(T, "j_r"),
        is_j_total_momentum = repeat([iszero(P.spin)], size(T)),
        is_calculated_with_mqdt = is_mqdt(T),
        underspecified_channel_contribution = get_neg(T),
    )
    return df
end


function exp_nui(T::MQDT.BasisArray)
    t = Vector{Float64}(undef, size(T))
    for (i, state) in enumerate(T.states)
        t[i] = exp_q(state.nu_list, state.coefficients)
    end
    return t
end

function std_nui(T::MQDT.BasisArray)
    t = Vector{Float64}(undef, size(T))
    for (i, state) in enumerate(T.states)
        t[i] = std_q(state.nu_list, state.coefficients)
    end
    return t
end

function is_mqdt(T::MQDT.BasisArray)
    t = Vector{Bool}(undef, size(T))
    for (i, state) in enumerate(T.states)
        t[i] = !isone(length(state.coefficients))
    end
    return t
end


function get_neg(T::MQDT.BasisArray)
    t = Vector{Float64}(undef, size(T))
    for (i, state) in enumerate(T.states)
        irrel = findall(iszero, state.model.core)
        t[i] = sum(state.coefficients[irrel] .^ 2)
    end
    return t
end

function exp_q(q::Vector, n::Vector)
    if allequal(q)
        return Float64(q[1])
    else
        m = n .^ 2
        M = sum(m)
        if M > 1
            m /= M
        end
        return LinearAlgebra.dot(q, m)
    end
end

function std_q(q::Vector, n::Vector)
    if allequal(q)
        return 0.0
    else
        m = n .^ 2
        M = sum(m)
        if M > 1
            m /= M
        end
        e1 = LinearAlgebra.dot(q, m)^2
        e2 = LinearAlgebra.dot(q .^ 2, m)
        if abs(e1 - e2) < 1e-11
            return 0.0
        else
            return sqrt(e2 - e1)
        end
    end
end

function get_n(T::MQDT.BasisArray, P::MQDT.Parameters)
    nu = MQDT.get_nu(T)
    l = round.(Int, calc_exp_qn(T, "l_r"))
    return get_n(nu, l, P.species)
end

function get_n(nu::Vector{Float64}, l::Vector{Int}, species::Symbol)
    i0 = findall(iszero, l)
    i1 = findall(iszero, l .- 1)
    i2 = findall(iszero, l .- 2)
    i3 = findall(iszero, l .- 3)
    j0 = findall(x->x<2, nu)
    nu[j0] .+= 1
    if occursin("Yb", String(species))
        nu[i0] .+= 4
        nu[i1] .+= 3
        nu[i2] .+= 2
        nu[i3] .+= 1
    else
        nu[i0] .+= 3
        nu[i1] .+= 2
        nu[i2] .+= 2
    end
    return ceil.(Int, nu)
end
