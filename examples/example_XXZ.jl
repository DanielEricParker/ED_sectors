######################################################################
#An example of a minimal version of the ED_sectors code

######################################################################

include("../main.jl")

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
println(abstract_H)

#make the Hamiltonian matrix
H = construct_matrix(basisFull, abstract_H)

#full ED
eigsys = full_ED(H)
println(typeof(eigsys.evs))
println(typeof(eigsys.O))

#show the eigenvalues
println(eigsys.evs[1:10])
