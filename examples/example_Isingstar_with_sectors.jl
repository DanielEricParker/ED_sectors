######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################

include("../main.jl")


############ XXZ* Hamiltonian Example with Sectors ###################

 function abstract_Ham_Ising_star(
 	L :: Int,
 	g_sigma :: Float64,
 	U_sigma :: Float64,
 	g_tau :: Float64,
 	U_tau :: Float64,
 	gamma :: Float64,
 	B_scale :: Float64)
	H = ABSTRACT_OP(L,"Ising_star",true)
	H += TERMS(-g_sigma,"ZIZ",0,2)
	H += TERMS(-U_sigma,"XIX",0,2)
	H += TERMS(-1.0, "X", 0, 2)

	H += TERMS(-B_scale, "X", 1,2)
	H += TERMS(-g_tau*B_scale, "ZIZ",1,2)
	H += TERMS(-g_tau*B_scale, "XIX",1,2)

	H += TERMS(-gamma, "XX", 0,2)

    return H
end

#set parameters
L = 12
a = 2

g_sigma = 1.0

U_sigma = 0.2
g_tau = 0.2
U_tau = 0.1
B_scale = 10.0
gamma = 0.0

println("Making abstract Hamiltonian...")
#make the abstract Hamiltonian from our parameters
abstract_hamiltonian = abstract_Ham_Ising_star(L,g_sigma,U_sigma,g_tau,U_tau,gamma,B_scale)

symmetries = Dict{String,Int64}("Z2A" => 1, "Z2B" => 1,"K"=>0) 
println("Making easy basis for symmetry sectors:\n\t",symmetries)
#make a basis with symmetries --- this can be a bit slow
basis = make_basis(L; unitCellSize=a, syms=symmetries)

println("Making full basis...")
#basis = make_basis(L)

println("Basis size: ", length(basis.conj_classes))

println("Making the Hamiltonian matrix...")
#make the Hamiltonian matrix
H = construct_matrix(basis,abstract_hamiltonian)


println("Finding smallest eigenvalues...")
#find the smallest 10 eigenvalues
evs = eigs(H; which=:SR, nev=30)

#show the eigenvalues
println(real(evs[1]))