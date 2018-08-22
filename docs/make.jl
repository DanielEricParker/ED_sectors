push!(LOAD_PATH,"/home/dan/Documents/ED_sectors/src/")

using Documenter
using ED_sectors

# makedocs(
# 	    # options
#     modules = [ED_sectors],
#     doctest = true
#     )

# deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
#     repo = "github.com/DanielEricParker/ED_sectors.git",
#     julia  = "1.0",
#     osname = "linux")


makedocs(
# 	    # options
    modules = [ED_sectors],
    format = :html,
#     doctest = true
	sitename = "ED_sectors",
	pages = [
		"index.md",
		"Abstract Operators" => "abstract_operators.md",
		"Basis" => "basis.md",
		"Matrix Constructors" => "matrix_constructors.md",
		"Measurement" => "measurement.md",
		"Dynamics" => "dynamics.md",
		"Utilities" => "utilities.md"
		],
	)


deploydocs(
    repo = "github.com/DanielEricParker/ED_sectors.git",
	target = "build",
    julia = "1.0",
    deps = nothing,
	make = nothing,
)