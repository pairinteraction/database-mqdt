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
include("tables.jl")
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
        "--n-min"
        help = "The minimal principal quantum number n for the states to be included in the database."
        default = 25
        arg_type = Int
        "--n-max"
        help = "The maximum principal quantum number n for the states to be included in the database."
        default = 90
        arg_type = Int
        "--directory"
        help = "The directory where the database will be saved"
        default = "database"
        arg_type = String
        "--overwrite"
        help = "Delete the species folder if it exists and create a new one"
        action = :store_true
        "--skip-high-l"
        help = "Include high angular momentum (l) states in the calculation"
        action = :store_true
    end

    return parse_args(s)
end


function main()
    # Parse command line arguments
    args = parse_commandline()
    n_min = args["n-min"]
    n_max = args["n-max"]
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
    @info "Parameters: n_min=$n_min, n_max=$n_max"
    start_time = time()

    # initialize Wigner symbol calculation
    if args["skip-high-l"]
        CGcoefficient.wigner_init_float(max(FMODEL_MAX_L[species], 5), "Jmax", 9)
    else
        CGcoefficient.wigner_init_float(n_max - 1, "Jmax", 9)
    end
    parameters = PARA_TABLE[species]

    @info "Calculating low ℓ MQDT states..."
    models = MODELS_TABLE[species]
    @timelog states = [eigenstates(n_min, n_max, M, parameters) for M in models]

    if args["skip-high-l"]
        @info "Skipping high ℓ states."
    else
        @info "Calculating high ℓ SQDT states..."
        l_max = n_max - 1
        l_start = FMODEL_MAX_L[species] + 1
        high_l_models = single_channel_models(species, l_start:l_max, parameters)
        @timelog high_l_states =
            [eigenstates(n_min, n_max, M, parameters) for M in high_l_models]
        states = vcat(states, high_l_states)
        models = vcat(models, high_l_models)
    end

    basis = basisarray(states, models)
    @info "Generated state table with $(length(basis.states)) states"

    @info "Converting states to database table..."
    states_df = basis_to_df(basis, parameters)

    @info "Calculating matrix elements..."
    @timelog row_col_value_dict = all_matrix_element(basis, parameters)

    @info "Converting matrix elements to database table..."
    matrix_elements_df_dict = Dict(k => rcv_to_df(v) for (k, v) in row_col_value_dict)

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
