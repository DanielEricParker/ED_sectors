################## MEASUREMENT.JL ##################################
#
#	This file contains functions used to measure observables
#	in the ground state or other single eigenstates of a Hamiltonian.
#
#	These methods can take full advantage of sparsity.
#
#
####################################################


"""

	k_eigvals(H,k)

Computes the `k` lowest-energy eigenvalues for the given Hamiltonian `H`. 

!!! note 
    
    This is a simple wrapper for Julia's `eigs` function, which uses the Lanczos algorithm. While this works with either sparse or full matrices for `H`, sparse matrices are much faster. The cost of the algorithm increases with `k` and should be only used when ``k \\ll \\operatorname{dim}(H)``. If you only need the ground state, then just use `k=1`.

"""
function k_eigvals(H :: AbstractMatrix, k :: Int)
		L = size(H)[1]
		k = min(k, div(L,2))#can't take too many with a small matrix, because Lanczos breaks
	return real(eigs(H; which=:SR, nev=k)[1])
end 

"""
	k_eigsys(H, k)

Computes the `k` lowest-energy eigenvalues and eigenvectors for the given Hamiltonian `H`. 

"""
function k_eigsys(H :: AbstractMatrix, k :: Int)
	L = size(H)[1]
	k = min(k, div(L,2))
	eig_sys =  eigs(H; which=:SR, nev=k)
	return (real(eig_sys[1]),eig_sys[2])
end 

#thought: for a single-site operator, it's probably way 
#better to compute the reduced density matrix
"""

	expectation(psi, O)
		expectation_rho(rho, O)

For a state `psi`, computes the expectation value ``\\left<\\psi|O|\\psi\\right>``  and operator `O`. This is a simple wrapper for `dot(psi, O*psi)` and works for any type of matrix `O`.

For a density matrix `rho`, computes the expectation value ``\\operatorname{Tr}[\\rho O]``.

"""
function expectation(
	psi :: Array{Complex{Float64}},
	Op :: AbstractMatrix
	)

	return dot(psi,Op*psi)

	#to-do
	#check if it's faster with BLAS or julia. probably depends on the type
	#<psi|Op|psi>
	# len = length(psi)
	# v = BLAS.gemv('N',Op,psi) #N for no transpose
	# corr = BLAS.dotc(len,psi,1,v,1)
	# return corr
end

function expectation(
	rho :: Array{Complex{Float64}},
	Op :: AbstractMatrix
	)
	return tr(Op*rho)
end



# #thought: for a single-site operator, it's probably way 
# #better to compute the reduced density matrix
# """
# Function to return a correlation (no time evolution)

# #Argument 
# * 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
# * 'SOP :: Array{Complex{Float64},2}': the operator to measure
# """
# function correlation1pt(
# 	psi :: Array{Complex{Float64}},
# 	op :: SOP,
# 	factor :: Float64,
# 	basis :: BASIS
# 	)
# 	len = length(psi)

# 	abstract_O = ABSTRACT_OP(basis.L,op)
# 	Op_mat = construct_matrix(basis, abstract_O)

# 	#<psi|Op|psi>
# 	v = Op_mat * psi #N for no transpose
# 	corr = BLAS.dotc(len,psi,1,v,1)

# 	return corr
# end



# """
# Function to compute n-point correlation functions for a state psi.

# """
# function correlation_npt(
# 	psi :: Array{Complex{Float64}},
# 	Ops :: Array{SOP},
# 	factor :: Float64,
# 	basis :: BASIS
# 	)

# end




#eventually put in functionality that works for any subset of the states,
#not just contiguous ones
"""

	reduce_density_matrix(phi, l)

Computes the reduced density matrix for the state `psi` by tracing out all but the first `l` spins.

!!! warning "Basis restriction"

    At present this only works when the basis has no symmetries. Support for reduced density matrices with symmetries is coming soon.

"""
function reduce_density_matrix(
	phi :: AbstractVector,#no type to support both sparse and normal vectors: fix
	l :: Int,
	basis :: BASIS
	)
	@assert length(phi) == 2^basis.L "The length of psi does not match the basis."
	@assert basis.q_numbers == Dict{String,Int}() " At present the reduced density matrix can only be computed for bases without symmetries. Support for reduced density matrices with symmetries is coming soon."

	new_dim = 2^l
	phi_reshape = reshape(phi,div(length(phi),new_dim),new_dim)

	println(typeof(phi_reshape))
	rho_reduced = [
			#come back and think if this should actually always be real
			real(dot(phi_reshape[:,i],phi_reshape[:,j]))
		for i in 1:new_dim, j in 1:new_dim]

	for i in 1:new_dim, j in 1:new_dim 
		if isnan(rho_reduced[i,j])

			println(phi_reshape[:,i])
			println(phi_reshape[:,j])
			println(dot(phi_reshape[:,i],phi_reshape[:,j]))
			error("NaN at (i,j) = (", i, ", ", j,")")
		end
	end
	return Hermitian(rho_reduced)
end


#untested
"""
	entanglement_entropy(psi, l; epsilon = 10e-15)

Computes the entanglement entropy ``S(l/L)`` for a state `psi`. Removes eigenvalues smaller than a numerical cutoff because they often give underflow errors. 

This uses `reduce_density_matrix`, so it only works without symmetries.

"""
function entanglement_entropy(
	psi :: AbstractVector, #no type to support both sparse and normal vectors: fix
	l :: Int,
	basis :: BASIS,
	epsilon = 10e-15#hardcoded tolerance for numerical error
	)

	rho = Matrix(reduce_density_matrix(psi,l,basis))

	println(size(rho))
	#println(rho)
	eigsys = EigSys(rho)

	#println(real(evs.values))

	S = Float64(0.0)
	tol = 10e-14#hardcoded tolerance for numerical error

	for p in real(eigsys.evs)
		if p > epsilon
			S += -p*log(p)
			#println(S)
		elseif p < -epsilon
			error("Negative eigenvalue in the reduced density matrix!")
		end
	end
	return S
end

"""

	measure_EE(psi, basis)

Computes the entanglement entropy at all cuts in a wavefunction `psi` in a basis `basis`. 
"""
function measure_EE(
	psi :: AbstractVector,
	basis :: BASIS)
	L = basis.L
	EE = Array{Float64,2}(undef,L+1,2)

	for k in 1:(floor(Int,L/2)+1)
		ee_k = entanglement_entropy(psi,k,basis;epsilon = 10e-15)
		EE[k,1] = k
		EE[k,2] = ee_k
		#this might overlap on the middle, but that's fine, it just overwrites the float
		EE[L+2-k,1] = L+2-k
		EE[L+2-k,2] = ee_k
	end

	return EE
end




# """
# Numerically estimates where a transition is by detecting a gapless spectrum with CFT scaling. Works for a 1 parameter family of Hamiltonians.
# #Argument 
# * 'L1 :: Int': the smaller size of Hamiltonian to try
# * 'L2 :: Int': the larger size of Hamiltonian to try
# * 'guess :: Float64': guess for the parameter where the transition occurs
# * 'Ham': 1 parameter family of Hamiltonians
# """
# function transition_finder(
# 	L1 :: Int,
# 	L2 :: Int,
# 	guess :: Float64,
# 	Ham
# 	)

# end