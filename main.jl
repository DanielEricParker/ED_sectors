######################################################################
#The main file for our project. Import this to get all the functionality

######################################################################

using Nullables
using Printf
using IterativeEigensolvers
using DelimitedFiles
using SparseArrays
using LinearAlgebra

#Kludgey way to get tests
testing = false
testing2 = false
testing3 = false
testing4 = false
testing5 = false

##################### Import everything #################


include("utilities.jl")
include("easy_basis.jl")
include("abstract_operators.jl")
include("matrix_constructors.jl")
include("measurement.jl")
include("dynamics.jl")