using MQDT
using CGcoefficient
using Parquet2
using ArgParse
import Base.Filesystem: mkpath
using DataFrames
using OrderedCollections
using Logging, LoggingExtras


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
        default = 20
        arg_type = Int
        "--n-max"
        help = "The maximum principal quantum number n for the states to be included in the database."
        default = 80
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
    @time states = [eigenstates(n_min, n_max, M, parameters) for M in models]

    if args["skip-high-l"]
        @info "Skipping high ℓ states."
    else
        @info "Calculating high ℓ SQDT states..."
        l_max = n_max - 1
        l_start = FMODEL_MAX_L[species] + 1
        high_l_models = single_channel_models(l_start:l_max, parameters)
        @time high_l_states =
            [eigenstates(n_min, n_max, M, parameters) for M in high_l_models]
        states = vcat(states, high_l_states)
        models = vcat(models, high_l_models)
    end

    basis = basisarray(states, models)
    state_table = state_data(basis, parameters)
    @info "Generated state table with $(nrow(state_table)) states"

    @info "Calculating matrix elements..."
    @time (d1, d2, dm, dd) = all_matrix_element(basis, parameters)

    @info "Converting matrix elements to database table..."
    m1 = matrix_data(d1)
    m2 = matrix_data(d2)
    mm = matrix_data(dm)
    md = matrix_data(dd)

    @info "Preparing database output..."
    db = databasearray(states, models)
    st = state_data(db, parameters)

    @info "Storing database tables as parquet files..."
    tables = OrderedDict(
        "states" => (data = st, desc = "States table"),
        "matrix_elements_d" => (data = m1, desc = "Dipole matrix elements"),
        "matrix_elements_q" => (data = m2, desc = "Quadrupole matrix elements"),
        "matrix_elements_mu" => (data = mm, desc = "Magnetic matrix elements"),
        "matrix_elements_q0" => (data = md, desc = "Diamagnetic matrix elements"),
    )
    for (name, table) in tables
        @info "$(table.desc) info" rows=nrow(table.data)
        @info describe(table.data)
        @time Parquet2.writefile("$(output_dir)/$(name).parquet", table.data)
    end

    elapsed_time = round(time() - start_time, digits = 2)
    @info "Database generation completed" elapsed_seconds=elapsed_time
    @info "Output saved to: $output_dir"

end


main()
