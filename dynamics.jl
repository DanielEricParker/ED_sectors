######################################################
#This file provides methods for dynamics with ED, inclduing
# - an "evolver" function for evolving states forward in time
# - a "correlation" function for measuring 1 and 2 pt fcns
# - a "timeseries" function for measuring a correlation
# 	with linear or log-spacing in time
# - a "thermal density matrix" function for computing Tr[rho_T O]
#These depend, of course, on full diagonalization.
#Eventually there should be a transparent interface   
#to do this with sectors all "under the hood".
#####################################################

###########Issue! I mixed up U and U^d
#Note: D  = U^d H U or H = U D U^d
#so exp(-i t H) = U exp(-it D) U^d
#incorrect functions are commentsed for now

"""
Struct for eigensystems, containing all eigenvectors and eigenvalues of matrix, the output of full ED.
"""
struct EigSys
	evs :: Array{Float64} #eigenvalues
	O :: Array{Complex{Float64},2} #eigenvectors


	"""
	Function to perform full diagonalization on the Hamiltonian.
	Returns an EigSys.
	"""
	function EigSys(H :: Matrix)
		eigensystem = eigfact(Matrix(H))
		new(real(eigensystem.values),eigensystem.vectors)
	end


	"""
	Alternatively, construct some other way.
	Returns an EigSys.
	"""
	function EigSys(
				evs :: Array{Float64}, #eigenvalues
				O :: Array{Complex{Float64},2} #eigenvectors
			)
		new(evs,O)
	end
end
Base.show(io::IO, eigsys :: EigSys) = print(io, "Evolver[] of size ", length(eigsys.evs))


"""
Function to return a correlation at a certain type.

#Arguments
* 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
* 'evolve :: Evolver': an evolver for our Hamiltonian
* 't :: Float64': the time to measure at
* 'Op :: OP': the operator to measure
* 'factor :: Float64': the scale on the operator
"""
#thought: for a single-site operator, it's probably way 
#better to compute the reduced density matrix
#for now, let's leverage the Hamiltonian constructor to 
#make a full operator, even though that's a bit dumb
#
#in particular, this makes no use of sparsity
function correlation1pt(
	psi :: Array{Complex{Float64},1},
	eigsys :: EigSys,
	t :: Float64,
	Op :: OP,
	factor :: Float64,
	basis :: Basis
	)

	abstract_O = ABSTRACT_OP(basis.L,Op)
	Op_mat = construct_matrix(basis, abstract_O)

	return correlation(psi,eigsys,t,Op_mat)
end

"""
Function to return a correlation at a certain type.

#Arguments
* 'psi :: Array{ComplexF64}': the (normalized) state vector to measure
* 'evolve :: Evolver': an evolver for our Hamiltonian
* 't :: Float64': the time to measure at
* 'Op :: Array{Complex{Float64},2}': the operator to measure
"""
#thought: for a single-site operator, it's probably way 
#better to compute the reduced density matrix

function correlation(
	psi :: Array{Complex{Float64}},
	eigsys :: EigSys,
	t :: Float64,
	Op :: AbstractMatrix
	)
	len = length(psi)
	

	expEigenvals1 = [exp(-im * t * ev) for ev in eigsys.evs]
	expEigenvals2 = [exp(im * t * ev) for ev in eigsys.evs]

	D1 = Diagonal(expEigenvals1)
	D2 = Diagonal(expEigenvals2)

	#<psi|Op(t)|psi>
	#	= <psi|exp(-itH) * Op * exp(itH)|psi>
	#	= <psi|U*exp(-itD)*U^d*Op*U*exp(itD)*U^d|psi> 

	#to make this maximally efficient, one should 
	#take advantage of the sparseness of the diagonal matrices
	#but let's not bother for now
	#
	#I'm not sure how to memory-manage this,
	#so for now I'm just over-writing v many times

	v = BLAS.gemv('C',eigsys.O,psi) #N for no transpose
	
	v = D2 * v #this is faster than -> v = BLAS.gemv('N',D2,v) 
	v = BLAS.gemv('N',eigsys.O,v) #C for conjugate transpose
	v = Op * v #faster in the case we have a sparse matrix
	#v = BLAS.gemv('N',Op, v)
	v = BLAS.gemv('C',eigsys.O,v) #N for no transpose
	v = D1 * v
	v = BLAS.gemv('N',eigsys.O,v) #C for conjugate transpose

	corr = BLAS.dotc(len,psi,1,v,1)

	return corr
end


"""
Function to return a correlation at a certain type.
* 'rho :: Array{Complex{Float64},2}': the (normalized) density matrix
* 'evolve :: Evolver': an evolver for our Hamiltonian
* 't :: Float64': the time to measure at
* 'Op :: Array{Complex{Float64},2}': the operator to measure

"""
function correlation1pt(
	rho :: Array{Complex{Float64},2},
	eigsys :: EigSys,
	t :: Float64,
	Op :: Array{Complex{Float64},2}
	)

	expEigenvals1 = [exp(-im * t * ev) for ev in eigsys.evs]
	expEigenvals2 = [exp(im * t * ev) for ev in eigsys.evs]

	D1 = Diagonal(expEigenvals1)
	D2 = Diagonal(expEigenvals2)

	# Tr[ rho U*exp(-itD)*U^d*Op*U*exp(itD)*U^d ]
	M = D2 * adjoint(eigsys.O) #D2 * U^d
	#M = BLAS.gemm('C','N',D2,evolve.O)
	M = BLAS.gemm('N','N',eigsys.O,M) # U *
	M = BLAS.gemm('N','N',Op,M)	# Op *
	M = BLAS.gemm('C','N',eigsys.O,M) #U^d *
	M = D1* M #exp(-itD) *
	#M = BLAS.gemm('C','N',D1,M)
	M = BLAS.gemm('N','N',eigsys.O,M) #U *
	M = BLAS.gemm('N','N',rho,M) #rho *

	return trace(M)

end	


"""
Function to return a correlation at a certain type.
* 'rho :: Array{Complex{Float64},2}': the (normalized) density matrix
* 'evolve :: Evolver': an evolver for our Hamiltonian
* 't :: Float64': the time to measure at
* 'Op :: Array{Complex{Float64},2}': the first operator
* 'Op :: Array{Complex{Float64},2}': the second operator
"""
#this is essentially the same as before
#perhaps I should just combine them
function correlation2pt(
	rho :: Array{Complex{Float64},2},
	eigsys :: EigSys,
	t :: Float64,
	Op1 :: Array{Complex{Float64},2},
	Op2:: Array{Complex{Float64},2}
	)


	expEigenvals1 = [exp(-im * t * ev) for ev in eigsys.evs]
	expEigenvals2 = [exp(im * t * ev) for ev in eigsys.evs]

	D1 = Diagonal(expEigenvals1)
	D2 = Diagonal(expEigenvals2)

	# Tr[ rho U*exp(-itD)*U^d*Op_1*U*exp(itD)*U^d Op_2]

	M = BLAS.gemm('C','N',eigsys.O, Op2)
	M = D2 * M #D2 * U
	#M = BLAS.gemm('C','N',D2,evolve.O)
	M = BLAS.gemm('N','N',eigsys.O,M) # U *
	M = BLAS.gemm('N','N',Op1,M)	# Op *
	M = BLAS.gemm('C','N',eigsys.O,M) #U^d *
	M = D1* M #exp(-itD) *
	#M = BLAS.gemm('C','N',D1,M)
	M = BLAS.gemm('N','N',eigsys.O,M) #U *
	M = BLAS.gemm('N','N',rho,M) #rho *

	return trace(M)
end

"""
Function to compute the thermal density matrix
"""
function thermal_density_matrix(
	eigsys :: EigSys,
	beta :: Float64
	)
	
	boltzmann = [exp(-beta * ev) for ev in eigsys.evs]
	Z = sum(boltzmann)
	D = (1/Z) * Diagonal(boltzmann)

	M = D * adjoint(eigsys.O)
	M = BLAS.gemm('N','N',eigsys.O,M)
	return M
end

#for infinite temeperature, i.e. beta = 0
function thermal_density_matrix(
	eigsys :: EigSys,
	beta :: Int
	)
	if(beta != 0)
		return thermal_density_matrix(eigsys,Float64(beta))
	end
	Z =length(eigsys.evs)
	a = 1/Z
	D = Matrix{ComplexF64}(Diagonal([a for i in 1:Z]))
	return D
end


"""
Function to compute a two-point function <psi|O(t)O(0)|psi>
at a list of times {t1,t2,...}. Since we're doing the same
operator at each time, we can optimize this a bit.
"""
function timeseries(
	psi :: Array{Complex{Float64}},
	eigsys :: EigSys,
	times :: Array{Float64},
	Op1 :: Array{Complex{Float64},2},
	Op2:: Array{Complex{Float64},2};
	verbose :: Bool = false
)
	len = length(psi)
	# <psi| Op_1(t) Op_2(0)|psi>
	# <psi| U D_m U^d Op_1 U D_p U^d Op_2|psi>
	# (<psi| U) * D_m * (U^d Op_1 U) * D_p * (U^d Op_2 |psi>)

	UdOp1U = BLAS.gemm('N','N',Op1,eigsys.O)
	UdOp1U = BLAS.gemm('C','N',eigsys.O,UdOp1U)

	UdOp2psi =  BLAS.gemv('N',Op2,psi) 
	UdOp2psi =  BLAS.gemv('C',eigsys.O,UdOp2psi) 

	Udpsi = BLAS.gemv('C',eigsys.O,psi)

	corr_t = Array{Float64,2}(uninitialized, length(times),3)
	for k in 1:length(times)
		t = times[k]
		if verbose
			println("Time $(k)/$(length(times))")
		end
		expEigenvals_m = [exp(-im * t * ev) for ev in eigsys.evs]
		expEigenvals_p = [exp(im * t * ev) for ev in eigsys.evs]

		D_m = Diagonal(expEigenvals_m)
		D_p = Diagonal(expEigenvals_p)

		v = D_p * UdOp2psi
		v = BLAS.gemv('N',UdOp1U,v)
		v = D_m * v

		cor = BLAS.dotc(len,Udpsi,1,v,1)
		corr_t[k,1] = t
		corr_t[k,2] = real(cor)
		corr_t[k,3] = imag(cor)
	end

	return corr_t
end

"""
Function to compute a two-point function Tr[rho O(t)O(0)]
at a list of times {t1,t2,...}. 
"""
function timeseries(
	rho :: Array{Complex{Float64},2},
	eigsys :: EigSys,
	times :: Array{Float64},
	Op1 :: Array{Complex{Float64},2},
	Op2:: Array{Complex{Float64},2};
	verbose :: Bool = false
)

	#we want 
	# Tr[ rho Op_1(t) Op_2(0)]
	# Tr[ Op_2 rho Op_1(t)]
	# Tr[ Op_2 rho    U exp(-it D) U^d  Op_1 U exp(it D) U^d]
	# Tr[ exp(itD) * (U^d Op2 rho U) * exp(-itD) * (U^d*Op1*U)  ]

	UdOp1U = BLAS.gemm('N','N',Op1,eigsys.O)
	UdOp1U = BLAS.gemm('C','N',eigsys.O,UdOp1U)

	UdOp2rhoU = BLAS.gemm('N','N',rho,eigsys.O)
	UdOp2rhoU = BLAS.gemm('N','N',Op2,UdOp2rhoU)
	UdOp2rhoU = BLAS.gemm('C','N',eigsys.O,UdOp2rhoU)

	corr_t = Array{Float64}(uninitialized, length(times), 3)
	for k in 1:length(times)
		t = times[k]
		if verbose
			println("Time $(k)/$(length(times))")
		end
		expEigenvals_m = [exp(-im * t * ev) for ev in eigsys.evs]
		expEigenvals_p = [exp(im * t * ev) for ev in eigsys.evs]

		D_m = Diagonal(expEigenvals_m)
		D_p = Diagonal(expEigenvals_p)

		M = D_p * UdOp2rhoU * D_m * UdOp1U
		cor = trace(M)#trace(M)

		corr_t[k,1] = t
		corr_t[k,2] = real(cor)
		corr_t[k,3] = imag(cor)
	end

	return corr_t

end


"""
Function to compute a two-point function <psi|O(t)O(0)|psi>
for an *eigenstate* |psi>
at a list of times {t1,t2,...}. Since we're doing the same
operator at each time, we can optimize this a bit.
"""
#turns out there should be a much better way to do this with many fewer
#matrix multiplications
function timeseries2(
	psi :: Array{Complex{Float64}}, ###must be an eigenstate
	E_psi :: Float64,
	eigsys :: EigSys,
	times :: Array{Float64},
	Op, #can be sparse
	verbose :: Bool = false
)
	len = length(psi)

	Oppsi = Op * psi

	O_sr = BLAS.gemv('C',eigsys.O,Oppsi)
	O_sr_abs = Array{Complex{Float64}}([abs2(x) for x in O_sr])


	corr_t = Array{Float64,2}(uninitialized, length(times),3)
	for k in 1:length(times)
		t = times[k]
		if verbose
			println("Time $(k)/$(length(times))")
		end
		exp_Es = exp(-im*t*E_psi)
		expEigenvals = [exp(im * t * ev) *exp_Es for ev in eigsys.evs]

		cor = BLAS.dotu(len,O_sr_abs,1,expEigenvals,1)
		corr_t[k,1] = t
		corr_t[k,2] = real(cor)
		corr_t[k,3] = imag(cor)
	end

	return corr_t
end