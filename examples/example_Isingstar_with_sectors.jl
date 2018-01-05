######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################

testing5 = false

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
	H = HAMILTONIAN("Ising_star",[])
    for i in 0:2:L-1
		H += term(-g_sigma, 			OP("Z", i), 	OP("Z", (i+2) % L))
		H += term(-U_sigma, 			OP("X", i), 	OP("X", (i+2) % L))
	 	H += term(-1.0,			 		OP("X", i))
	 	H += term(-1.0*B_scale, 		OP("X", i+1))
		H += term(-g_tau*B_scale, 		OP("Z", i+1), 	OP("Z", (i+3) % L))
		H += term(-U_tau*B_scale, 		OP("X", i+1), 	OP("X", (i+3) % L))
		H += term(-gamma, 				OP("X", i), 	OP("X", (i+1) % L))
    end
    return H
end

#set parameters
L = 12
a = 2

g_sigma = 1.0
U_sigma = 0.0
g_tau = 0.0
U_tau = 0.0
B_scale = 1.0
gamma = 0.0

# U_sigma = 0.2
# g_tau = 0.2
# U_tau = 0.1
# B_scale = 10.0
# gamma = 0.0

println("Making abstract Hamiltonian...")
#make the abstract Hamiltonian from our parameters
abstract_hamiltonian = abstract_Ham_Ising_star(L,g_sigma,U_sigma,g_tau,U_tau,gamma,B_scale)

validity_symmetries = Dict{String,Int64}()#Dict("SzA" => 3)
symmetries = Dict{String,Int64}("Z2A" => -1, "Z2B" => -1,"K"=>1) 
println("Making easy basis for symmetry sectors:\n\t", validity_symmetries,"\n\t",symmetries)
#make a basis with symmetries --- this can be a bit slow
SzTrbasis2 = make_easy_basis(L,a,validity_symmetries,symmetries)

# @time SzTrbasis2 = make_universal_basis(L,a,validity_symmetries,symmetries)
# @time SzTrbasis2 = make_universal_basis(L,a,validity_symmetries,symmetries)

# x = UInt64(134)
# println("x:", bin(x))

# y = flip_spins_B([x],8)
# println(bin(y[1]))
# orbit1 = get(apply_Z2B(8,1,ORBIT(1,x,Dict(x => 1.0))))

# orb_fcn = make_easy_orbit_function(8,2,Dict("Z2B" => 1))

# orbit2 = get(orb_fcn(x))

# println("rep1:", bin(orbit1.representative))
# for (k,v) in orbit1.elements
# 	println(bin(k), ", ", v)
# end

# println("rep2:", bin(orbit2.representative))
# for (k,v) in orbit2.elements
# 	println(bin(k), ", ", v)
# end

# println("Making full basis...")
# SzTrbasis2 = make_full_basis(L)

println("Basis size: ", length(SzTrbasis2.conj_classes))

println("Making the Hamiltonian matrix...")
#make the Hamiltonian matrix
H2 = make_Hamiltonian(L,SzTrbasis2,abstract_hamiltonian)

println("Finding smallest eigenvalues...")
#find the smallest 10 eigenvalues
evs2 = eigs(H2;which=:SR, nev=30)

#show the eigenvalues
println(real(evs2[1]))
