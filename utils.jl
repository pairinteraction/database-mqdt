import SparseArrays: spzeros
import MQDT

function all_matrix_element(B::BasisArray, parameters::MQDT.Parameters)
    """Calculate all relevant matrix elements for a given basis array B.

    This means dipole, quadrupole, magnetic, and diamagnetic matrix elements.
    """
    k_angular_max = 2  # 2 for now, since we dont calculate octupole matrix elements

    n_states = length(B.states)
    dipole = spzeros(n_states, n_states)
    quadrupole = spzeros(n_states, n_states)
    magnetic = spzeros(n_states, n_states)
    diamagnetic = spzeros(n_states, n_states)

    list_of_qns =
        [(ids, state.nu, state.lr, state.f) for (ids, state) in enumerate(B.states)]
    qns_sorted_by_lr = sort(list_of_qns, by = x -> (minimum(x[3]), minimum(x[2]), x[1]))
    ids_sorted = [x[1] for x in qns_sorted_by_lr]
    # ids_sorted = [i for i in 1:n_states]

    for (i1, id1) in enumerate(ids_sorted)
        b1 = B.states[id1]
        f1 = b1.f

        for id2 in ids_sorted[i1:end]
            b2 = B.states[id2]
            f2 = b2.f

            # TODO
            # if all(n > all_n_up_to for n in [n1, n2]) && abs(n1 - n2) > max_delta_n
            #     continue
            # end

            # dipole
            value = MQDT.multipole_moment(1, b1, b2)
            if !iszero(value)
                dipole[id1, id2] = value
                if id2 != id1
                    dipole[id2, id1] = (-1)^(f2 - f1) * value
                end
            end

            # quadrupole
            value = MQDT.multipole_moment(2, b1, b2)
            if !iszero(value)
                quadrupole[id1, id2] = value
                if id2 != id1
                    quadrupole[id2, id1] = (-1)^(f2 - f1) * value
                end
            end

            # magnetic
            value = MQDT.magnetic_dipole_moment(parameters.dipole, parameters.mass, parameters.spin, b1, b2)
            if !iszero(value)
                magnetic[id1, id2] = value
                if id2 != id1
                    magnetic[id2, id1] = (-1)^(f2 - f1) * value
                end
            end

            # diamagnetic
            value = MQDT.special_quadrupole_moment(b1, b2)
            if !iszero(value)
                diamagnetic[id1, id2] = value
                if id2 != id1
                    diamagnetic[id2, id1] = (-1)^(f2 - f1) * value
                end
            end

        end
    end
    return (dipole, quadrupole, magnetic, diamagnetic)
end
