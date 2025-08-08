using MQDT
using CGcoefficient
using Parquet2
using ArgParse
import Base.Filesystem: mkpath

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
        default = 30
        arg_type = Int
        "--n-max"
        help = "The maximum principal quantum number n for the states to be included in the database."
        default = 40
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

# Parse command line arguments
args = parse_commandline()
n_min = args["n-min"]
n_max = args["n-max"]
directory = args["directory"]
overwrite = args["overwrite"]
species = Symbol(args["species"])

# initialize Wigner symbols
CGcoefficient.wigner_init_float(n_max, "Jmax", 9) # initialize Wigner symbol caluclation

parameters = PARA_TABLE[species]

# calculate low \ell MQDT states
low_l_models = MODELS_TABLE[species]
low_l_states = [eigenstates(n_min, n_max, M, parameters) for M in low_l_models]

# calculate high \ell SQDT states
l_max = n_max - 1
high_l_models = single_channel_models(5:l_max, parameters)
high_l_states = [eigenstates(n_min, n_max, M, parameters) for M in high_l_models]

# generate state table
basis = basisarray(vcat(low_l_states, high_l_states), vcat(low_l_models, high_l_models))
state_table = state_data(basis, parameters)

# calculate matrix elements
@time d1 = matrix_element(1, basis) # dipole
@time d2 = matrix_element(2, basis) # quadrupole
@time dm = matrix_element(parameters, basis) # Zeeman
@time dd = matrix_element(basis) # diamagnetic

# convert matrix elements to database table
m1 = matrix_data(d1)
m2 = matrix_data(d2)
mm = matrix_data(dm)
md = matrix_data(dd)

# prepare database output
db = databasearray(vcat(low_l_states, high_l_states), vcat(low_l_models, high_l_models))
ST = state_data(db, parameters)

# store full matrix (as opposed to upper triangle)
M1 = tri_to_full(m1, ST)
M2 = tri_to_full(m2, ST)
MM = tri_to_full(mm, ST)
MD = tri_to_full(md, ST)


# store tables as parquet files
output_dir = "$(directory)/$(species)_mqdt_$(version)"
if overwrite
    rm(output_dir, recursive = true, force = true)
elseif isdir(output_dir)
    error("Output directory already exists: $output_dir. Use --overwrite to overwrite.")
end
mkpath(output_dir)
Parquet2.writefile("$(output_dir)/states.parquet", ST)
Parquet2.writefile("$(output_dir)/matrix_elements_d.parquet", M1)
Parquet2.writefile("$(output_dir)/matrix_elements_q.parquet", M2)
Parquet2.writefile("$(output_dir)/matrix_elements_mu.parquet", MM)
Parquet2.writefile("$(output_dir)/matrix_elements_q0.parquet", MD)
