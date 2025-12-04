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

    states_indexed = [(ids - 1 + START_ID, state) for (ids, state) in enumerate(B.states)]
    states_sorted =
        sort(states_indexed, by = x -> (minimum(x[2].lr), minimum(x[2].nu), x[1]))

    for (i1, (id1, b1)) in enumerate(states_sorted)
        for (id2, b2) in states_sorted[i1:end]

            # Skip if all contributions of the two states are far apart in angular momentum
            if minimum(b2.lr) - maximum(b1.lr) > k_angular_max
                continue
            end

            # Skip if all contributions of the two states are far apart in n and None of them is low-n
            if all(abs(nu1-nu2) >= 11 for nu1 in b1.nu for nu2 in b2.nu) &&
               all(nu1 > 25 for nu1 in b1.nu) &&
               all(nu2 > 25 for nu2 in b2.nu)
                continue
            end

            m = MQDT.multipole_moments(b1, b2, parameters)
            # multipole_moments returns the matrix elements in the following order
            # electric dipole, electric quadrupole, diamagnetic, magnetic
            table_keys = [
                "matrix_elements_d",
                "matrix_elements_q",
                "matrix_elements_q0",
                "matrix_elements_mu",
            ]
            prefactor_transposed = (-1)^(b2.f - b1.f)

            for (i, key) in enumerate(table_keys)
                if m[i] != 0
                    push!(row_col_value[key], (id1, id2, m[i]))
                    if id1 != id2
                        push!(row_col_value[key], (id2, id1, m[i] * prefactor_transposed))
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


function get_n(T::MQDT.BasisArray)
    # TODO for now just return round(nu)
    # later calculate radial overlap with different sqdt states and take the corresponding n
    nu = MQDT.get_nu(T)
    return round.(Int, nu)
end


function basis_to_df(T::MQDT.BasisArray, P::MQDT.Parameters)
    df = DataFrame(
        id = collect(START_ID:(size(T)-1+START_ID)),
        energy = MQDT.get_e(T, P) / 219474.6313632,  # convert 1/cm to atomic units
        parity = MQDT.get_p(T),
        n = get_n(T),
        nu = MQDT.get_nu(T),
        f = MQDT.get_f(T),
        exp_nui = MQDT.exp_nui(T),
        exp_l = calc_exp_qn(T, "l_tot"),
        exp_j = calc_exp_qn(T, "j_tot"),
        exp_s = calc_exp_qn(T, "s_tot"),
        exp_l_ryd = calc_exp_qn(T, "l_r"),
        exp_j_ryd = calc_exp_qn(T, "j_r"),
        std_nui = MQDT.std_nui(T),
        std_l = calc_std_qn(T, "l_tot"),
        std_j = calc_std_qn(T, "j_tot"),
        std_s = calc_std_qn(T, "s_tot"),
        std_l_ryd = calc_std_qn(T, "l_r"),
        std_j_ryd = calc_std_qn(T, "j_r"),
        is_j_total_momentum = MQDT.is_J(T, P),
        is_calculated_with_mqdt = MQDT.is_mqdt(T),
        underspecified_channel_contribution = MQDT.get_neg(T),
    )
    return df
end
