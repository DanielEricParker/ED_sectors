######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################

include("../main.jl")


############ XXZ* Hamiltonian Example with Sectors ###################

 function make_XXZ_star_operators_new(L :: Int, Delta :: Float64, alpha :: Float64, g_tau :: Float64, u_tau :: Float64, B_scale :: Float64)
	H = HAMILTONIAN("XXZ_star",[])
    for i in 0:2:L-1
		H += term(0.5, 				OP("+", i), 	OP("-", (i+2) % L))
		H += term(0.5, 				OP("-", i), 	OP("+", (i+2) % L))
		H += term(Delta, 			OP("Z", i), 	OP("Z", (i+2) % L))
	 	H += term(-1.0*B_scale, 	OP("X",(i+1) % L))
	 	H += term(-g_tau*B_scale, 	OP("Z",(i+1) % L),OP("Z",(i+3) % L))
	 	H += term(-u_tau*B_scale, 	OP("X",(i+1) % L),OP("X",(i+3) % L))
	 	H += term(alpha, 			OP("Z",i),		OP("X",(i+1) % L),	OP("Z", (i+2)%L))
    end
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


validity_symmetries = Dict("SzA" => 3)
symmetries = Dict{String,Int64}("K" => 3, "Z2B" => 1)
println("Making basis for symmetry sectors:\n\t", validity_symmetries,"\n\t",symmetries)
#make a basis with symmetries --- this can be a bit slow
SzTrbasis2 = make_universal_basis(L,a,validity_symmetries,symmetries)


println("Basis size: ", length(SzTrbasis2.conj_classes))

println("Making the Hamiltonian matrix...")
#make the Hamiltonian matrix
H2 = make_Hamiltonian(L,SzTrbasis2,abstract_hamiltonian)

println("Finding smallest eigenvalues...")
#find the smallest 10 eigenvalues
evs2 = eigs(H2;which=:SR, nev=10)

#show the eigenvalues
println(real(evs2[1]))