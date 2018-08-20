include("../main.jl")

using Test

@testset "Running all tests..." begin
	include("test_abstract_operators.jl")
	include("tests_easy_basis.jl")
	include("tests_matrix_constructors.jl")
end