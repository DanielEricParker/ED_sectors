######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################

include("../main.jl")


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
H = make_Hamiltonian(L, basisFull, abstract_H)

#full ED
evs = eigfact(Matrix(H))
println(typeof(evs.values))
println(typeof(evs.vectors))

#show the eigenvalues
println(evs.values[1:10])
