######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################


using BenchmarkTools, Compat

include("../main.jl")

# @btime [exp(2*x) for x in 1:1000000]

################### XXZ Hamiltonian example ##########################


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
H = Matrix(make_Hamiltonian(L, basisFull, abstract_H))

#full ED
evs = eigfact(Matrix(H))
println(typeof(evs.values))
println(typeof(evs.vectors))

# #show the eigenvalues
# println(evs.values[1:10])

#make the evolver
println("Testing the evolver...")
evolve = evolver(H)
println(evolve)

#test the eigenvalues and how to access them
n = 3
vec = evolve.O[:,n]	#get nth eigenvector
w = evolve.evs[n]	#get nth eigenvalue
v2 = H * vec
v3 = w * vec
println(norm(v2-v3)) #should be zero


#test the 1pt correlation fcn 
println("Testing the 1pt correlation fcn")
psi = evolve.O[:,25] #the ground state, because why not
t = 10.0
Op = OP("X",3)
factor = 1.0

println("Finding 1pt function from abstract Op")
corr1pt = correlation1pt(psi, evolve, t, Op, factor, basisFull)
println(corr1pt)

println("Finding 1pt function from matrix Op")
abstract_O = HAMILTONIAN("O", [term(factor,OP("Z",3),OP("Z",4))])
Op_mat = Matrix(make_Hamiltonian(L, basisFull, abstract_O))

corr = correlation(psi, evolve, t, Op_mat)
println(corr)

println("Finding 1pt function from matrix Op and density matrix rho")
rho = kron(transpose(psi),psi)
corr_rho = correlation1pt(rho, evolve, t, Op_mat)
println(corr_rho)


ts = [Float64(t) for t in 1:20]

println("Finding 2pt function for rho")
corr_2pt_rho = correlation2pt(rho,evolve,t,Op_mat,Op_mat)
println(corr_2pt_rho)

println("Finding 2pt timeseries for psi")
timesseries_rho = timeseries(psi,evolve,[t],Op_mat,Op_mat)
println(timesseries_rho)


println("Finding 2pt timeseries for rho")
timesseries_rho = timeseries(rho,evolve,[t],Op_mat,Op_mat)
println(timesseries_rho)



@btime [correlation2pt(rho,evolve,t1,Op_mat,Op_mat) for t1 in ts]
@btime timeseries(rho,evolve,ts,Op_mat,Op_mat) #obviously this is way faster
#it's about a factor 3 for 10 times
@btime timeseries(psi,evolve,ts,Op_mat,Op_mat) #and this is waaaay faster yet