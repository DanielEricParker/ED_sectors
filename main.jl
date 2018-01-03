######################################################################
#The main file for our project. Import this to get all the functionality

######################################################################

using Nullables
using IterativeEigensolvers

#Kludgey way to get tests
testing = false
testing2 = false
testing3 = false
testing4 = false

##################### Import everything #################


include("utilities.jl")
include("basis.jl")
include("abstract_Hamiltonian.jl")