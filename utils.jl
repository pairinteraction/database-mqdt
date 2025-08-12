import MQDT

function all_matrix_element(B::BasisArray, parameters::MQDT.Parameters)
    """Calculate all relevant matrix elements for a given basis array B.

    This means dipole, quadrupole, magnetic, and diamagnetic matrix elements.
    """
    k_angular_max = 2  # 2 for now, since we dont calculate octupole matrix elements

    row_col_value = Dict(
        "dipole" => Tuple{Int64,Int64,Float64}[],
        "quadrupole" => Tuple{Int64,Int64,Float64}[],
        "magnetic" => Tuple{Int64,Int64,Float64}[],
        "diamagnetic" => Tuple{Int64,Int64,Float64}[],
    )

    indexed_states = [(ids, state) for (ids, state) in enumerate(B.states)]
    sorted_states =
        sort(indexed_states, by = x -> (minimum(x[2].lr), minimum(x[2].nu), x[1]))

    for (i1, (id1, b1)) in enumerate(sorted_states)
        # TODO filter to include only states with similar l (via k_angular_max)

        for (id2, b2) in sorted_states[i1:end]
            # TODO
            # if all(n > all_n_up_to for n in [n1, n2]) && abs(n1 - n2) > max_delta_n
            #     continue
            # end

            value_dipole = MQDT.multipole_moment(1, b1, b2)
            if value_dipole != 0
                push!(row_col_value["dipole"], (id1, id2, value_dipole))
            end

            value_quadrupole = MQDT.multipole_moment(2, b1, b2)
            if value_quadrupole != 0
                push!(row_col_value["quadrupole"], (id1, id2, value_quadrupole))
            end

            value_magnetic = MQDT.magnetic_dipole_moment(
                parameters.dipole,
                parameters.mass,
                parameters.spin,
                b1,
                b2,
            )
            if value_magnetic != 0
                push!(row_col_value["magnetic"], (id1, id2, value_magnetic))
            end

            value_diamagnetic = MQDT.special_quadrupole_moment(b1, b2)
            if value_diamagnetic != 0
                push!(row_col_value["diamagnetic"], (id1, id2, value_diamagnetic))
            end

            if id1 != id2
                prefactor = (-1)^(b2.f - b1.f)
                if value_dipole != 0
                    push!(row_col_value["dipole"], (id2, id1, prefactor * value_dipole))
                end
                if value_quadrupole != 0
                    push!(
                        row_col_value["quadrupole"],
                        (id2, id1, prefactor * value_quadrupole),
                    )
                end
                if value_magnetic != 0
                    push!(row_col_value["magnetic"], (id2, id1, prefactor * value_magnetic))
                end
                if value_diamagnetic != 0
                    push!(
                        row_col_value["diamagnetic"],
                        (id2, id1, prefactor * value_diamagnetic),
                    )
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
