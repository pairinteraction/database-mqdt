# Tool for generating the MQDT database tables in the cloud

Database of states and matrix elements calculated with multichannel quantum defect theory.
Database tables are available through [GitHub Releases](https://github.com/pairinteraction/database-mqdt/releases).

## Generate new tables locally
To generate the database tables locally, first make sure you also downloaded MQDT.jl as submodule (e.g. run `git submodule update --init`).

Then, initialize the julia environment and link against the MQDT.jl package
```bash
julia --project=. -e 'import Pkg; Pkg.develop(path="MQDT.jl"); Pkg.instantiate()'
```

Finally, you can run the script to generate the tables via
```bash
julia --project=. generate_database.jl Yb174
```

## Generate a new release
To generate a new release with all the important tables, simply create and push a new annotated tag with a tag name of the form `v*.*` .
This will run the `generate_database.yml` workflow, where first all tables are created (this happens for all commits, not only tags),
and then in addition uploads the zipped versions of the tables to a new release with name `v*.*` .
The release is created in draft mode, so you can double-check, that all database tables are included and optionally add a release text.
Once you are happy with the release draft, don't forget to publish the release.
