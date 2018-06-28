####################################################
#functions for measurements
#should eventually include:
# -1pt, 2pt, npt functions (sparse)
# -arbitary observables (dense)
# -reduced density matrix
# -entanglement entropy
# -transition finder based on eigenvalue crossings (perhaps put in utilities instead?)
###################s#################################


"""
Function to compute the first k eigenvalues.
"""
function k_eigvals(H :: AbstractMatrix, k :: Int)
	L = size(H)[1]
	k = min(k, div(L,2))#can't take too many with a small matrix, because Lanczos breaks
	return real(eigs(H; which=:SR, nev=k)[1])
end 

"""
Function to compute the first k eigenvalues and eigenvectors.
"""
function k_eigsys(H :: AbstractMatrix, k :: Int)
	L = size(H)[1]
	k = min(k, div(L,2))
	return eigs(H; which=:SR, nev=k)
end 

"""
Function to perform full ED. Actually just a wrapper for the EigSys constructor.
"""
function full_ED(H :: AbstractMatrix)
	return EigSys(Matrix(H)) #convert to dense, because really there's no point in this for a sparse matrix. it's the same speed I believe.
end

"""
Function to make an observable with nice syntax.

#Argument 
* 'basis :: BASIS': the basis we're working with
* 'op :: String': the name for the operator we want to measure, e.g. "ZZ"
* 'site :: Int': starting site for the operator. 
"""
function observable(
	basis :: Basis,
	op :: String,
	site :: Int
	)
	abstract_op  = ABSTRACT_OP(basis.L,OP(op,site))
	op_mat = construct_matrix(basis, abstract_op)
	return op_mat
end




#thought: for a single-site operator, it's probably way 
#better to compute the reduced density matrix
"""
Function to return a correlation with no time evolution

#Argument 
* 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
* 'Op :: Array{Complex{Float64},2}': the operator to measure
"""
function correlation(
	psi :: Array{Complex{Float64}},
	Op :: Matrix
	)
	len = length(psi)

	#<psi|Op|psi>
	v = BLAS.gemv('N',Op,psi) #N for no transpose
	corr = BLAS.dotc(len,psi,1,v,1)

	return corr
end

#thought: for a single-site operator, it's probably way 
#better to compute the reduced density matrix
"""
Function to return a correlation (no time evolution)

#Argument 
* 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
* 'Op :: Array{Complex{Float64},2}': the operator to measure
"""
function correlation1pt(
	psi :: Array{Complex{Float64}},
	op :: OP,
	factor :: Float64,
	basis :: Basis
	)
	len = length(psi)

	abstract_O = ABSTRACT_OP(basis.L,op)
	Op_mat = construct_matrix(basis, abstract_O)

	#<psi|Op|psi>
	v = Op_mat * psi #N for no transpose
	corr = BLAS.dotc(len,psi,1,v,1)

	return corr
end



"""
Function to compute n-point correlation functions for a state psi.

"""
function correlation_npt(
	psi :: Array{Complex{Float64}},
	Ops :: Array{OP},
	factor :: Float64,
	basis :: Basis
	)

end



"""
Function to compute n-point correlation functions for a density matrix rho.

"""
function correlation_npt(
	rho :: Array{Complex{Float64},2},
	Ops :: Array{OP},
	factor :: Float64,
	basis :: Basis
	)

end



#eventually put in functionality that works for any subset of the states,
#not just contiguous ones
"""
Computes the reduced density matrix for a state psi.
#Argument 
* 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
* 'l :: Int': the number of sites we want to have left
"""
function reduce_density_matrix(
	phi,#no type to support both sparse and normal vectors: fix
	l :: Int
	)
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
Computes the entanglement entropy S(l/L) for a state psi. 
"""
function entanglement_entropy(
	psi, #no type to support both sparse and normal vectors: fix
	l :: Int;
	epsilon = 10e-15#hardcoded tolerance for numerical error
	)

	rho = Matrix(reduce_density_matrix(psi,l))

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

"""Computes EE for cuts at every site for a wavefunction.
#Argument 
* 'L :: Int': the number of sites
* 'psi :: Array{Complex{Float64},1}': the wavefunction
"""
function measure_EE(
	L :: Int,
	psi :: Array{Complex{Float64},1})
	EE = Array{Float64,2}(undef,L+1,2)

	for k in 1:(floor(Int,L/2)+1)
		ee_k = entanglement_entropy(psi,k;epsilon = 10e-15)
		EE[k,1] = k
		EE[k,2] = ee_k
		#this might overlap on the middle, but that's fine, it just overwrites the float
		EE[L+2-k,1] = L+2-k
		EE[L+2-k,2] = ee_k
	end

	return EE
end




"""
Numerically estimates where a transition is by detecting a gapless spectrum with CFT scaling. Works for a 1 parameter family of Hamiltonians.
#Argument 
* 'L1 :: Int': the smaller size of Hamiltonian to try
* 'L2 :: Int': the larger size of Hamiltonian to try
* 'guess :: Float64': guess for the parameter where the transition occurs
* 'Ham': 1 parameter family of Hamiltonians
"""
function transition_finder(
	L1 :: Int,
	L2 :: Int,
	guess :: Float64,
	Ham
	)

end

#verified against mathematica code
#this is easy to read rather than fast, but it shouldn't matter much
"""
Computes the R-statistic for a spectrum. Must be in a single symmetry sector. For integrable systems, R ~ 0.386. For GOE level-statistics, R ~0.5295. 
"""
function r_statistic(
	eigenvalues :: Array{Float64}
	)
	#must be sorted
	eigenvalues = sort(eigenvalues)
	rs = Array{Float64}(undef, length(eigenvalues)-2)

	for i in 1:length(rs)
		delta_na = eigenvalues[i+1]-eigenvalues[i]
		delta_nb = eigenvalues[i+2]-eigenvalues[i+1]
		r_min = min(delta_nb,delta_na)
		r_max = max(delta_nb,delta_na)
		r = (abs(r_min - r_max) < 10e-14) ? 0.0 : r_min/r_max
		#just in case, check we're in the bounds
		if r < -0.001 || r > 1.001
			error("out of bounds", i," ", r_min," ", r_max," ", r)
		end
		rs[i] = r
	end
	return mean(rs)
end