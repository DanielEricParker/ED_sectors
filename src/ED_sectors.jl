# __precompile__()

module ED_sectors

	using Base
	using Nullables
	using Printf
	using IterativeSolvers
	using DelimitedFiles
	using SparseArrays
	using LinearAlgebra
	using Arpack


	include("utilities.jl")
	export export_data, make_filename

	include("abstract_operators.jl")
	#export datatypes
	export TERM, ABSTRACT_OP
	#no functions to export!

	include("basis.jl")
	#export datatypes
	export BASIS
	#no functions to export! they're all internal
	#this should be revisited so users can add their own symmetries

	include("matrix_constructors.jl")
	#export datatypes --- none needed
	#export functions
	export Operator

	include("measurement.jl")
	#no datatypes to export!
	#export functions
	export k_eigvals, k_eigsys, full_ED, observable, correlation, correlation1pt, reduce_density_matrix, entanglement_entropy, measure_EE, r_statistic

	include("dynamics.jl")
	#export datatypes
	export EigSys
	#export functions
	export evolve_state, correlation1pt, correlation, correlation1pt, correlation2pt, thermal_density_matrix, timeseries,timeseries2

end	