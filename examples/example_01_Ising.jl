######################################################################
#example_01_Ising
#An example of a minimal version of the ED_sectors code

######################################################################

#include the ED_Sectors codebase
include("../src/ED_sectors.jl")



#Define a function to create a symbolic Ising Hamiltonian
"""

	symbolic_Ising_Hamiltonian(L, g)

A function to generate the (symbolic) Ising Hamiltonian with periodic boundary conditions.

# Arguments
- `L::Int`: the number of sites in the spin chain.
- `g::Float64`: parameter for the Ising model. g = 0 is a pure ferromagnet, g -> infinity is a pure paramagnet, and g = 1 is the transition.
	
"""
function symbolic_Ising_Hamiltonian(L :: Int, g :: Float64)
    H = ABSTRACT_OP(L, "Ising", true) #(number of sites, name, periodic boundaries?)
    H += TERMS(g,"X")	#sum_i g*sigma_x^i
    H += TERMS(1.0,"ZZ") #sum_i sigma_z^i sigma_z^{i+1}
    return H
end


#A wrapper function for all our code
#Julia is JIT compiled, so we get better performance if
#we wrap our code in a function so it compiles first
function main()
	println("Setting parameters...")
	L = 8
	g = 0.8
	println("L = $(L), g= $(0.8)")

	println("\nMaking a basis for a spin chain with $(L) sites, with no symmetries for now")
	basisFull = Basis(L)
	println(basisFull)

	println("\nMaking the symbolic Hamiltonian...")
	symbolic_H = symbolic_Ising_Hamiltonian(L,g)
	println(symbolic_H)


	println("\nMaking the Hamiltonian matrix from the symbolic one...")
	H = construct_matrix(basisFull, symbolic_H)

	println("\nH is a matrix of size $(size(H)).")

	println("\nPerforming full diagonalization...")
	eigsys = full_ED(H)


	println("\nThe lowest 10 eigenvalues are:")
	println(eigsys.evs[1:10])
end

main()
