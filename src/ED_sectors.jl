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
	using Statistics

	include("utilities.jl")
	export export_data, make_filename

	include("abstract_operators.jl")
	#export datatypes
	export TERM, ABSTRACT_OP,TERM_INTERNAL,consolidate_coefficients
	#no functions to export!

	include("basis.jl")
	#export datatypes
	export BASIS
	#no functions to export! they're all internal
	#this should be revisited so users can add their own symmetries

	include("pauli_strings.jl")
	export VERT, PAULI, showVerbose, commutator, shift_left, shift_right,mult, PauliOp,adj,trim,string_length

	include("matrix_constructors.jl")
	#export datatypes --- none needed
	#export functions
	export Operator,PauliString,construct_matrix_full,PauliOp,construct_matrix_full_new,construct_matrix_full_new_Pauli,to_Pauli_strings

	include("measurement.jl")
	#no datatypes to export!
	#export functions
	export k_eigvals, k_eigsys, expectationValue, reduce_density_matrix, entanglement_entropy, measure_EE

	include("full_spectrum_measurement.jl")
	#export datatypes
	export EIGSYS
	export thermal_density_matrix, r_statistic


	include("dynamics.jl")

	#export functions
	export evolve_state, expectation_time, autocorrelation
	#, correlation1pt, correlation, correlation1pt, correlation2pt, thermal_density_matrix, timeseries,timeseries2

end	