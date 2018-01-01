######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################

testing3 = false


include("main.jl")


###############Make a Hamiltonian constructor  #######################


function abstract_XXZ_Hamiltonian(L :: Int, Delta :: Float64)
    #L -- length of the spin chain
    #Delta -- anisotropy parameter
    H = HAMILTONIAN("XXZ", [])
    for i in 0:L-1
        H += term(0.5, OP("+",i), OP("-",(i+1)%L))
        H += term(0.5, OP("-",i), OP("+",(i+1)%L))
        H += term(Delta*1.0, OP("Z",i), OP("Z",(i+1)%L))
    end

    return H
end

#set parameters
L = 8
Delta = 0.8

#make the basis
basisFull = make_full_basis(L)

#make the abstract Hamiltonian
abstract_H = abstract_XXZ_Hamiltonian(L,Delta)

#make the Hamiltonian matrix
H = make_Hamiltonian(L, basisFull, abstract_H)

#diagonalize
evs = eigfact(Matrix(H))

println(evs.values[1:10])



################ test universal basis function

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


L = 16

Delta = 0.8
alpha = 0.1
g_tau  = 0.3
u_tau = 0.1
B_scale = 1.0

@time SzTrbasis1 = make_Sz_Tr_basis(L,4,1)

#println("Basis size:", length(SzTrbasis1.conj_classes))
#println(collect(Set(values(SzTrbasis1.get_conj_class))))

@time SzTrbasis2 = make_universal_basis(L,2,Dict("SzA" => 4,"K" => 1))

#println("Basis size:", length(SzTrbasis2.conj_classes))
#println(collect(Set(values(SzTrbasis2.get_conj_class))))


abstract_hamiltonian = make_XXZ_star_operators_new(L,Delta,alpha,g_tau,u_tau,B_scale)


H1 = make_Hamiltonian(L,SzTrbasis1,abstract_hamiltonian)
evs1 = eigfact(Matrix(H1))
println(evs1.values[1:10])

H2 = make_Hamiltonian(L,SzTrbasis2,abstract_hamiltonian)
evs2 = eigfact(Matrix(H2))
println(evs2.values[1:10])