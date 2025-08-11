using MQDT
using CGcoefficient
using Parquet2
using ArgParse
import Base.Filesystem: mkpath
using DataFrames
using OrderedCollections


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


    println("Starting database generation for $species with version $version")
    println("Parameters: n_min=$n_min, n_max=$n_max")
    start_time = time()

    CGcoefficient.wigner_init_float(n_max, "Jmax", 9) # initialize Wigner symbol calculation
    parameters = PARA_TABLE[species]

    println("Calculating low ℓ MQDT states...")
    low_l_models = MODELS_TABLE[species]
    @time low_l_states = [eigenstates(n_min, n_max, M, parameters) for M in low_l_models]

    println("Calculating high ℓ SQDT states...")
    l_max = n_max - 1
    high_l_models = single_channel_models(5:l_max, parameters)
    @time high_l_states = [eigenstates(n_min, n_max, M, parameters) for M in high_l_models]

    basis = basisarray(vcat(low_l_states, high_l_states), vcat(low_l_models, high_l_models))
    state_table = state_data(basis, parameters)
    println("Generated state table with $(nrow(state_table)) states")

    println("Calculating matrix elements...")
    @time d1 = matrix_element(1, basis) # dipole
    @time d2 = matrix_element(2, basis) # quadrupole
    @time dm = matrix_element(parameters, basis) # Zeeman
    @time dd = matrix_element(basis) # diamagnetic

    println("Converting matrix elements to database table...")
    m1 = matrix_data(d1)
    m2 = matrix_data(d2)
    mm = matrix_data(dm)
    md = matrix_data(dd)

    println("Preparing database output...")
    db = databasearray(vcat(low_l_states, high_l_states), vcat(low_l_models, high_l_models))
    st = state_data(db, parameters)

    println("Storing database tables as parquet files...")
    tables = OrderedDict(
        "states" => (data = st, desc = "States table"),
        "matrix_elements_d" => (data = m1, desc = "Dipole matrix elements"),
        "matrix_elements_q" => (data = m2, desc = "Quadrupole matrix elements"),
        "matrix_elements_mu" => (data = mm, desc = "Magnetic matrix elements"),
        "matrix_elements_q0" => (data = md, desc = "Diamagnetic matrix elements"),
    )
    for (name, table) in tables
        println("\n$(table.desc) info:")
        println("Number of rows: $(nrow(table.data))")
        println(describe(table.data))
        @time Parquet2.writefile("$(output_dir)/$(name).parquet", table.data)
    end

    elapsed_time = round(time() - start_time, digits = 2)
    println("Database generation completed after in total $elapsed_time seconds")
    println("Output saved to: $output_dir")

end

main()
