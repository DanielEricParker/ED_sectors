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
Function to return a correlation with no time evolution

#Argument 
* 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
* 'Op :: Array{Complex{Float64},2}': the operator to measure
"""
#thought: for a single-site operator, it's probably way 
#better to compute the reduced density matrix

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


"""
Function to return a correlation (no time evolution)

#Argument 
* 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
* 'Op :: Array{Complex{Float64},2}': the operator to measure
"""
#thought: for a single-site operator, it's probably way 
#better to compute the reduced density matrix

function correlation1pt(
	psi :: Array{Complex{Float64}},
	Op :: OP,
	factor :: Float64,
	basis :: Basis
	)
	#	this is dumb. perhaps Basis should contain L?
	L = Int(round(log(2,length(psi))))
	len = length(psi)

	abstract_O = HAMILTONIAN("O", [term(factor,Op)])
	Op_mat = make_Hamiltonian(L, basis, abstract_O)


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




"""
Computes the reduced density matrix for a state psi.
#Argument 
* 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
* 'l :: Int': the number of sites we want to have left
"""
#eventually put in functionality that works for any subset of the states,
#not just contiguous ones
function reduce_density_matrix(
	phi :: Array{Complex{Float64}},
	l :: Int
	)


end



"""
Computes the entanglement entropy S(l/L) for a state psi. 
"""
#untested
function entanglement_entropy(
	psi :: Array{Complex{Float64}},
	l :: Int
	)

	rho = reduce_density_matrix(phi,l)
	evs = eigfact(rho)

	S = Float64(0.0)
	tol = 10e-16#hardcoded tolerance for numerical error

	for p in evs.values
		if abs(p) > tol
			S += -1*p*log(p)
		elseif p < -tol
			error("Negative eigenvalue in the reduced density matrix!")
		end
	end
	return S
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


"""
Computes the R-statistic for a spectrum. Must be in a single symmetry sector. For integrable systems, R ~ 0.386. For GOE level-statistics, R ~0.5295. 
"""
#verified against mathematica code
#this is easy to read rather than fast, but it shouldn't matter much
function r_statistic(
	eigenvalues :: Array{Float64}
	)
	#must be sorted
	eigenvalues = sort(eigenvalues)
	rs = Array{Float64}(uninitialized, length(eigenvalues)-2)

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