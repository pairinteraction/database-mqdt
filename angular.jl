using PythonCall
using MQDT

function calc_exp_qn(basisstate::BasisState, qn::String)
    state = get_state(basisstate)
    value = state.calc_exp_qn(qn)
    return pyconvert(Float64, value)
end

function calc_exp_qn(basisarray::BasisArray, qn::String)
    exp_list = Vector{Float64}(undef, size(basisarray))
    for (i, state) in enumerate(basisarray.states)
        exp_list[i] = calc_exp_qn(state, qn)
    end
    return exp_list
end

function calc_std_qn(basisstate::BasisState, qn::String)
    state = get_state(basisstate)
    value = state.calc_std_qn(qn)
    return pyconvert(Float64, value)
end

function calc_std_qn(basisarray::BasisArray, qn::String)
    std_list = Vector{Float64}(undef, size(basisarray))
    for (i, state) in enumerate(basisarray.states)
        std_list[i] = calc_std_qn(state, qn)
    end
    return std_list
end

function get_state(basisstate::BasisState)
    rydstate = pyimport("rydstate")
    channels = basisstate.channels.i
    kets = Vector{Any}(undef, size(channels))
    for i in eachindex(channels)
        qn = channels[i]
        kets[i] = get_ket(qn, basisstate.species)
    end
    coeff = basisstate.coefficients[findall(basisstate.model.core)]
    state = rydstate.angular.AngularState(coeff, kets; warn_if_not_normalized = false)
    return state
end

function get_ket(qn::lsQuantumNumbers, species::Symbol)
    rydstate = pyimport("rydstate")
    ket = rydstate.angular.AngularKetLS(;
        s_c = qn.sc,
        s_tot = qn.S,
        l_c = qn.lc,
        l_r = qn.lr,
        l_tot = qn.L,
        j_tot = qn.J,
        f_tot = qn.F,
        species = String(species),
    )
    return ket
end

function get_ket(qn::jjQuantumNumbers, species::Symbol)
    rydstate = pyimport("rydstate")
    ket = rydstate.angular.AngularKetJJ(;
        s_c = qn.sc,
        l_c = qn.lc,
        l_r = qn.lr,
        j_c = qn.Jc,
        j_r = qn.Jr,
        j_tot = qn.J,
        f_tot = qn.F,
        species = String(species),
    )
    return ket
end

function get_ket(qn::fjQuantumNumbers, species::Symbol)
    rydstate = pyimport("rydstate")
    ket = rydstate.angular.AngularKetFJ(;
        s_c = qn.sc,
        l_c = qn.lc,
        l_r = qn.lr,
        j_c = qn.Jc,
        f_c = qn.Fc,
        j_r = qn.Jr,
        f_tot = qn.F,
        species = String(species),
    )
    return ket
end
