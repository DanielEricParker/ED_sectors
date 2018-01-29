######################################################################
#example_02_XXZ.jl
#An example of the XXZ Hamiltonian, using sectors

######################################################################

include("../main.jl")

#Define a function to create the (symbolic) XXZ Hamiltonian
"""

	symbolic_XXZ_Hamiltonian(L, Delta)

A function to generate the (symbolic) XXZ Hamiltonian with periodic boundary conditions.

# Arguments
- `L::Int`: the number of sites in the spin chain.
- `Delta::Float64`: the anisotropy parameter for the XXZ model.

"""
function symbolic_XXZ_Hamiltonian(L :: Int, Delta :: Float64)
    #ABSTRACT_OP(number of sites, name, periodic boundary conditions?)
    H = ABSTRACT_OP(L, "XXZ", true)
    #H = sum_i sigma_{+}^i sigma_{-}^{i+1}
    H += TERMS(0.5, "+-")
    #H = sum_i sigma_{-}^i sigma_{+}^{i+1}
    H += TERMS(0.5, "-+")
    #H = sum_i Delta*sigma_z^i sigma_z^{i+1}
    H += TERMS(Delta*1.0, "ZZ")
    return H
end


#A wrapper function for all our code
#Julia is JIT compiled, so we get better performance if
#we wrap our code in a function so it compiles first
function main()
	println("Setting parameters...")
	L = 12
	g = 0.8
	println("L = $(L), g= $(0.8)")

	#Select the Sz and K symmetry sectors.
	#Sz is the number of spin ups, i.e. Sz = Mz - L/2
	#K is the translation sector
	#this sector contains the ground state
	symmetries = Dict{String,Int64}("Sz" => 6, "K" => 0)
	println("\nSymmetries are:")
	println(symmetries)

	println("\nMaking a basis for a spin chain with $(L) sites, with no symmetries for now")
	basis = make_basis(L;syms=symmetries)
	println("Basis is a basis of size $(basis.L).")


	println("\nMaking the symbolic Hamiltonian...")
	symbolic_H = symbolic_XXZ_Hamiltonian(L,g)
	println(symbolic_H)


	println("\nMaking the Hamiltonian matrix from the symbolic one...")
	H = construct_matrix(basis, symbolic_H)

	println("\nH is a matrix of size $(size(H)).")

	println("\nPerform Laczos diagonalization for the first 10 eiganvalues")
	eigs = k_eigvals(H,10)


	println("\nThe lowest 10 eigenvalues are:")
	println(eigs)
end

main()