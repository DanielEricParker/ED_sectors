######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################

include("../main.jl")


############ XXZ* Hamiltonian Example with Sectors ###################

 function make_XXZ_star_operators_new(L :: Int, Delta :: Float64, alpha :: Float64, g_tau :: Float64, u_tau :: Float64, B_scale :: Float64)
	H = ABSTRACT_OP(L,"XXZ_star",true)
	H += TERMS(0.5,"+I-",0,2)
	H += TERMS(0.5,"-I+",0,2)
	H += TERMS(Delta,"ZIZ",0,2)

	H += TERMS(-1.0*B_scale, "X",1,2)
	H += TERMS(-g_tau*B_scale, "ZZ", 1,2)
	H += TERMS(-u_tau*B_scale, "XX", 1,2)

	H += TERMS(alpha, "ZXZ", 0,2)

    return H
end

#set parameters
L = 12
a = 2

Delta = 0.8
alpha = 0.1
g_tau  = 0.3
u_tau = 0.1
B_scale = 1.0

println("Making abstract Hamiltonian...")
#make the abstract Hamiltonian from our parameters
abstract_hamiltonian = make_XXZ_star_operators_new(L,Delta,alpha,g_tau,u_tau,B_scale)

symmetries = Dict{String,Int64}("SzA" => 3,"K" => 3, "Z2B" => 1)
println("Making basis for symmetry sectors: \n\t", symmetries)
#make a basis with symmetries --- this can be a bit slow
basis = make_basis(L; unitCellSize = a, syms=symmetries)


println("Basis size: ", length(basis.conj_classes))

println("Making the Hamiltonian matrix...")
#make the Hamiltonian matrix
H = construct_matrix(basis,abstract_hamiltonian)

println("Finding smallest eigenvalues...")
#find the smallest 10 eigenvalues
evs = k_eigvals(H,10)

#show the eigenvalues
println(evs)