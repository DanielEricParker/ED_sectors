######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################


using BenchmarkTools#, Compat

include("../main.jl")

# @btime [exp(2*x) for x in 1:1000000]

################### XXZ Hamiltonian example ##########################


function abstract_XXZ_Hamiltonian(L :: Int, Delta :: Float64)
    #L -- length of the spin chain
    #Delta -- anisotropy parameter
    H = ABSTRACT_OP(L, "XXZ", true)
    H += TERMS(0.5, "+-")
    H += TERMS(0.5, "-+")
    H += TERMS(Delta*1.0, "ZZ")

    return H
end

#set parameters
L = 8
Delta = 0.8

#make the basis
basisFull = make_basis(L)

#make the abstract Hamiltonian
abstract_H = abstract_XXZ_Hamiltonian(L,Delta)

#make the Hamiltonian matrix
H = Matrix(construct_matrix(basisFull, abstract_H))

# #full ED
# evs = eigfact(Matrix(H))
# println(typeof(evs.values))
# println(typeof(evs.vectors))

# #show the eigenvalues
# println(evs.values[1:10])

#make the eigsysr
println("Testing the eigsys...")
eigsys = EigSys(H)
println(eigsys)

#test the eigenvalues and how to access them
n = 3
vec = eigsys.O[:,n]	#get nth eigenvector
w = eigsys.evs[n]	#get nth eigenvalue
v2 = H * vec
v3 = w * vec
println(norm(v2-v3)) #should be zero


#test the 1pt correlation fcn 
#println("Testing the 1pt correlation fcn")
psi = eigsys.O[:,25] #the 25th state, because why not
E_psi = eigsys.evs[25]
t = 10.0
ts = [Float64(t) for t in 1:20]
Op = OP("X",3)
factor = 1.0

println("Finding 1pt function from abstract Op")
corr1pt = correlation1pt(psi, eigsys, t, Op, factor, basisFull)
println(corr1pt)

#println("Finding 1pt function from matrix Op")
abstract_O = ABSTRACT_OP(L,"O",false,[TERM(factor,[OP("Z",3),OP("Z",4)])])
Op_mat_sp = construct_matrix(basisFull, abstract_O)
Op_mat = Matrix(Op_mat_sp)

corr = correlation(psi, eigsys, t, Op_mat)
println(corr)

println("Finding 1pt function from matrix Op and density matrix rho")
rho = kron(transpose(psi),psi)
corr_rho = correlation1pt(rho, eigsys, t, Op_mat)
println(corr_rho)


ts = [Float64(t) for t in 1:20]

println("Finding 2pt function for rho")
corr_2pt_rho = correlation2pt(rho,eigsys,t,Op_mat,Op_mat)
println(corr_2pt_rho)

println("Finding 2pt timeseries for psi")
timesseries_rho = timeseries(psi,eigsys,[t],Op_mat,Op_mat)
println(timesseries_rho)


println("Finding 2pt timeseries for rho")
timesseries_rho = timeseries(rho,eigsys,[t],Op_mat,Op_mat)
println(timesseries_rho)



@btime [correlation2pt(rho,eigsys,t1,Op_mat,Op_mat) for t1 in ts]
@btime timeseries(rho,eigsys,ts,Op_mat,Op_mat) #obviously this is way faster
#it's about a factor 3 for 10 times
@btime timeseries(psi,eigsys,ts,Op_mat,Op_mat) #and this is waaaay faster yet



#####second generation functions, which should be way more efficient

#old version
println("Finding 2pt timeseries for psi")
timesseries_psi = timeseries(psi,eigsys,ts,Op_mat,Op_mat)
println(timesseries_psi)
#@btime timeseries(psi,eigsys,ts,Op_mat,Op_mat)

#old version
println("Finding 2pt timeseries for psi --- faster?")
timesseries_psi_2 = timeseries2(psi,E_psi,eigsys,ts,Op_mat_sp)
println(timesseries_psi_2)
#@btime timeseries2(psi,E_psi,eigsys,ts,Op_mat_sp)

println("Are they approximately equal? ", isapprox(timesseries_psi,timesseries_psi_2))