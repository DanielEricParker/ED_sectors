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
		"Overview" => "index.md",
		"Tutorial" => "tutorial.md",
		"Full API" =>
			[
			"Abstract Operators" => "abstract_operators.md",
			"Making Bases" => "basis.md",
			"Operators" => "matrix_constructors.md",
			"Ground State Measurements" => "measurement.md",
			"Full Spectrum Measurements" => "full_spec.md",
			"Dynamics" => "dynamics.md",
			"Utilities" => "utilities.md",
			],
		"Internals" => "internals.md"],
	)


deploydocs(
    repo = "github.com/DanielEricParker/ED_sectors.git",
	target = "build",
    julia = "1.0",
    deps = nothing,
	make = nothing,
)