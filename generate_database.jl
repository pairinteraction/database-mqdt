using MQDT
using CGcoefficient
using Parquet2
using ArgParse
import Base.Filesystem: mkpath
using DataFrames
using OrderedCollections
using Logging, LoggingExtras
using LRUCache


include("utils.jl")
HIGH_L = 10

version = "v1.1"

function parse_commandline()
    s = ArgParseSettings(
        description = "Generate a database, containing energies and matrix elements, for a given species.",
        epilog = "Example:\n  julia --project=. generate_database.jl Yb174",
    )

    @add_arg_table! s begin
        "species"
        help = "The species to generate the database for"
        required = true
        arg_type = String
        "--n-min-sqdt"
        help = "The minimal principal quantum number n for the SQDT models to be included in the database."
        default = 25
        arg_type = Int
        "--n-max-high-l"
        help = "The maximal principal quantum number n for high angular momentum states (l > $HIGH_L) to be included in the database."
        default = Inf
        arg_type = Float64
        "--n-max"
        help = "The maximum principal quantum number n for the states to be included in the database."
        default = 110
        arg_type = Int
        "--directory"
        help = "The directory where the database will be saved"
        default = "database"
        arg_type = String
        "--overwrite"
        help = "Delete the species folder if it exists and create a new one"
        action = :store_true
    end

    return parse_args(s)
end


function main()
    # Parse command line arguments
    args = parse_commandline()
    n_max = args["n-max"]
    n_min_sqdt = args["n-min-sqdt"]
    n_max_high_l = args["n-max-high-l"]
    directory = args["directory"]
    overwrite = args["overwrite"]
    species = Symbol(args["species"])

    # Setup output directory
    output_dir = "$(directory)/$(species)_mqdt_$(version)"
    if overwrite
        rm(output_dir, recursive = true, force = true)
    elseif isdir(output_dir)
        error("Output directory already exists: $output_dir. Use --overwrite to overwrite.")
    end
    mkpath(output_dir)

    # Setup logging to both file and console
    console_logger = ConsoleLogger(stdout)
    file_logger = SimpleLogger(open("$(output_dir)/$(species).log", "w"))
    combined_logger = TeeLogger(console_logger, file_logger)
    global_logger(combined_logger)

    @info "Starting database generation for $species with version $version"
    @info "Parameters: $args"
    start_time = time()

    # initialize Wigner symbol calculation
    CGcoefficient.wigner_init_float(n_max - 1, "Jmax", 9)

    parameters = MQDT.get_species_parameters(species)
    all_models = Vector{MQDT.fModel}()

    s_r = 1 / 2
    j_c = 1 / 2
    i_c = parameters.spin
    for l_r = 0:(n_max-1)
        for j_r = abs(l_r-s_r):1:(l_r+s_r)
            for f_c = abs(j_c-i_c):1:(j_c+i_c)
                for f_tot = abs(f_c-j_r):1:(f_c+j_r)
                    models = MQDT.get_fmodels(species, l_r, j_r, f_c, f_tot, parameters)
                    for model in models
                        if !any(m -> m.name == model.name, all_models)
                            push!(all_models, model)
                        end
                    end
                end
            end
        end
    end

    @info "Calculating MQDT states..."
    states = Vector{MQDT.EigenStates}(undef, length(all_models))
    for (i, M) in enumerate(all_models)
        _n_min = NaN
        _n_max = n_max
        if startswith(M.name, "SQDT")
            _n_min = n_min_sqdt
            if parameters.spin > 0
                l_ryd = MQDT.get_lr(M)[1]
                if l_ryd >= n_max_high_l
                    _n_max = 0
                elseif l_ryd > HIGH_L
                    _n_max = min(n_max_high_l, n_max)
                end
            end
        end
        @info "$(M.name)"
        states[i] = MQDT.eigenstates(_n_min, _n_max, M, parameters)
        if length(states[i].n) > 0
            @info "  found nu_min=$(minimum(states[i].n)), nu_max=$(maximum(states[i].n)), total states=$(length(states[i].n))"
        else
            @info "  found no states"
        end
    end

    basis = MQDT.basisarray(states, all_models)
    @info "Generated state table with $(length(basis.states)) states"

    @info "Converting states to database table..."
    @timelog states_df = basis_to_df(basis, parameters)

    @info "Calculating matrix elements..."
    @timelog row_col_value_dict = all_matrix_element(basis, parameters)

    @info "Converting matrix elements to database table..."
    @timelog matrix_elements_df_dict =
        Dict(k => rcv_to_df(v) for (k, v) in row_col_value_dict)

    @info "Storing database tables as parquet files..."
    tables = merge(Dict("states" => states_df), matrix_elements_df_dict)
    for (name, table) in tables
        @info "$(name) info" rows=nrow(table)
        @info describe(table)
        @timelog Parquet2.writefile("$(output_dir)/$(name).parquet", table)
    end

    elapsed_time = round(time() - start_time, digits = 2)
    @info "Database generation completed" elapsed_seconds=elapsed_time
    @info "Output saved to: $output_dir"

    @info "Cache information MQDT.lru_radial_moment:"
    @info cache_info(MQDT.lru_radial_moment)
    @info "Cache information MQDT.lru_get_rydberg_state:"
    @info cache_info(MQDT.lru_get_rydberg_state)
    @info "Cache information MQDT.lru_angular_matrix:"
    @info cache_info(MQDT.lru_angular_matrix)
    @info "Cache information MQDT.lru_magnetic_matrix:"
    @info cache_info(MQDT.lru_magnetic_matrix)

end


main()
