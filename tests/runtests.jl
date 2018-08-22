#!/usr/bin/env julia

#Start Test Script
push!(LOAD_PATH,"../src/")

using ED_sectors
using Test

# Run tests

@testset "Running all tests..." begin
	include("test_abstract_operators.jl")
	include("tests_basis.jl")
	include("tests_matrix_constructors.jl")
end