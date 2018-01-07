######################################################
#This file provides methods for dynamics with ED, inclduing
# - an "evolver" function for evolving states forward in time
# - a "correlation" function for measuring 1 and 2 pt fcns
# - a "timeseries" function for measuring a correlation
# 	with linear or log-spacing in time
# - a "thermal density matrix" function for computing Tr[rho_T O]
#These depende, of course, on full diagonalization.
#Eventually there should be a transparent interface   
#to do this with sectors all "under the hood".
#####################################################




struct Evolver
	evs :: Array{Float64} #eigenvalues
	O :: Array{Complex{Float64},2} #eigenvectors
end
Base.show(io::IO, evolve :: Evolver) = print(io, "Evolver[] of size ", length(evolve.evs))


"""
Function to perform full diagonalization on the Hamiltonian.
Returns an Evolver (which is really just the eigenvalues and vectors).
"""
function evolver(H :: Matrix)
	eigensystem = eigfact(Matrix(H))
	return Evolver(real(eigensystem.values),eigensystem.vectors)
end



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
	evolve :: Evolver,
	t :: Float64,
	Op :: OP,
	factor :: Float64,
	basis :: Basis
	)
	#this is dumb. perhaps Basis should contain L?
	L = Int(round(log(2,length(psi))))
	#println("L: ", L)

	abstract_O = HAMILTONIAN("O", [term(factor,Op)])
	Op_mat = Matrix(make_Hamiltonian(L, basis, abstract_O))

	return correlation(psi,evolve,t,Op_mat)
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
	evolve :: Evolver,
	t :: Float64,
	Op :: Matrix
	)
	len = length(psi)
	

	expEigenvals1 = [exp(-im * t * ev) for ev in evolve.evs]
	expEigenvals2 = [exp(im * t * ev) for ev in evolve.evs]

	D1 = Diagonal(expEigenvals1)
	D2 = Diagonal(expEigenvals2)

	#<psi|Op(t)|psi>
	#	= <psi|exp(-itH) * Op * exp(itH)|psi>
	#	= <psi|U^d*exp(-itD)*U*Op*U^d*exp(itD)*U|psi> 

	#to make this maximally efficient, one should 
	#take advantage of the sparseness of the diagonal matrices
	#but let's not bother for now
	#
	#I'm not sure how to memory-manage this,
	#so for now I'm just over-writing v many times

	v = BLAS.gemv('N',evolve.O,psi) #N for no transpose
	
	v = D2 * v
	#v = BLAS.gemv('N',D2,v) #this is slower
	v = BLAS.gemv('C',evolve.O,v) #C for conjugate transpose
	v = BLAS.gemv('N',Op, v)
	v = BLAS.gemv('N',evolve.O,v) #N for no transpose
	#v = BLAS.gemv('N',D1,v)#this is slower
	v = D1 * v
	v = BLAS.gemv('C',evolve.O,v) #C for conjugate transpose

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
	evolve :: Evolver,
	t :: Float64,
	Op :: Array{Complex{Float64},2}
	)

	expEigenvals1 = [exp(-im * t * ev) for ev in evolve.evs]
	expEigenvals2 = [exp(im * t * ev) for ev in evolve.evs]

	D1 = Diagonal(expEigenvals1)
	D2 = Diagonal(expEigenvals2)

	# Tr[ rho U^d*exp(-itD)*U*Op*U^d*exp(itD)*U ]
	M = D2 * evolve.O #D2 * U
	#M = BLAS.gemm('C','N',D2,evolve.O)
	M = BLAS.gemm('C','N',evolve.O,M) # U^d *
	M = BLAS.gemm('N','N',Op,M)	# Op *
	M = BLAS.gemm('N','N',evolve.O,M) #U *
	M = D1* M #exp(-itD) *
	#M = BLAS.gemm('C','N',D1,M)
	M = BLAS.gemm('C','N',evolve.O,M) #U^d *
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
	evolve :: Evolver,
	t :: Float64,
	Op1 :: Array{Complex{Float64},2},
	Op2:: Array{Complex{Float64},2}
	)


	expEigenvals1 = [exp(-im * t * ev) for ev in evolve.evs]
	expEigenvals2 = [exp(im * t * ev) for ev in evolve.evs]

	D1 = Diagonal(expEigenvals1)
	D2 = Diagonal(expEigenvals2)

	# Tr[ rho U^d*exp(-itD)*U*Op_1*U^d*exp(itD)*U Op_2]

	M = BLAS.gemm('N','N',evolve.O, Op2)
	M = D2 * M #D2 * U
	#M = BLAS.gemm('C','N',D2,evolve.O)
	M = BLAS.gemm('C','N',evolve.O,M) # U^d *
	M = BLAS.gemm('N','N',Op1,M)	# Op *
	M = BLAS.gemm('N','N',evolve.O,M) #U *
	M = D1* M #exp(-itD) *
	#M = BLAS.gemm('C','N',D1,M)
	M = BLAS.gemm('C','N',evolve.O,M) #U^d *
	M = BLAS.gemm('N','N',rho,M) #rho *

	return trace(M)
end

"""
Function to compute the thermal density matrix
"""
####needs to be tested
function thermal_density_matrix(
	evolve :: Evolver,
	beta :: Float64
	)
	
	boltzmann = [exp(-beta * ev) for ev in evolve.evs]
	Z = sum(boltzmann)
	D = (1/Z) * Diagonal(boltzmann)

	M = D * evolve.O
	M = BLAS.gemm('C',evolve.O,M)
	return M
end

"""
Function to compute a two-point function <psi|O(t)O(0)|psi>
at a list of times {t1,t2,...}. Since we're doing the same
operator at each time, we can optimize this a bit.
"""
function timeseries(
	psi :: Array{Complex{Float64}},
	evolve :: Evolver,
	times :: Array{Float64},
	Op1 :: Array{Complex{Float64},2},
	Op2:: Array{Complex{Float64},2}
)
	len = length(psi)
	# <psi| Op_1(t) Op_2(0)|psi>
	# <psi| U^d D_m U Op_1 U^d D_p U Op_2|psi>
	# (<psi| U^d) * D_m * (U Op_1 U^d) * D_p * (U Op_2 |psi>)

	UOp1Ud = BLAS.gemm('N','C',Op1,evolve.O)
	UOp1Ud = BLAS.gemm('N','N',evolve.O,UOp1Ud)

	UOp2psi =  BLAS.gemv('N',Op2,psi) 
	UOp2psi =  BLAS.gemv('N',evolve.O,UOp2psi) 

	Upsi = BLAS.gemv('N',evolve.O,psi)

	corr_t = Array{Complex{Float64}}(uninitialized, length(times))
	for k in 1:length(times)
		t = times[k]
		expEigenvals_m = [exp(-im * t * ev) for ev in evolve.evs]
		expEigenvals_p = [exp(im * t * ev) for ev in evolve.evs]

		D_m = Diagonal(expEigenvals_m)
		D_p = Diagonal(expEigenvals_p)

		v = D_p * UOp2psi
		v = BLAS.gemv('N',UOp1Ud,v)
		v = D_m * v

		corr_t[k] = BLAS.dotc(len,Upsi,1,v,1)
	end

	return corr_t
end

"""
Function to compute a two-point function Tr[rho O(t)O(0)]
at a list of times {t1,t2,...}. 
"""
function timeseries(
	rho :: Array{Complex{Float64},2},
	evolve :: Evolver,
	times :: Array{Float64},
	Op1 :: Array{Complex{Float64},2},
	Op2:: Array{Complex{Float64},2}
)

	#we want 
	# Tr[ rho Op_1(t) Op_2(0)]
	# Tr[ Op_2 rho Op_1(t)]
	# Tr[ Op_2 rho    U^d exp(-it D) U  Op_1 U^d exp(it D) U]
	# Tr[ exp(itD) * (U Op2 rho U^d) * exp(-itD) * (U*Op1*U^d)  ]

	UOp1Ud = BLAS.gemm('N','C',Op1,evolve.O)
	UOp1Ud = BLAS.gemm('N','N',evolve.O,UOp1Ud)

	UOp2rhoUd = BLAS.gemm('N','C',rho,evolve.O)
	UOp2rhoUd = BLAS.gemm('N','N',Op2,UOp2rhoUd)
	UOp2rhoUd = BLAS.gemm('N','N',evolve.O,UOp2rhoUd)

	corr_t = Array{Complex{Float64}}(uninitialized, length(times))
	for k in 1:length(times)
		t = times[k]
		expEigenvals_m = [exp(-im * t * ev) for ev in evolve.evs]
		expEigenvals_p = [exp(im * t * ev) for ev in evolve.evs]

		D_m = Diagonal(expEigenvals_m)
		D_p = Diagonal(expEigenvals_p)

		M = D_p * UOp2rhoUd * D_m * UOp1Ud
		corr_t[k] = trace(M)#trace(M)
	end

	return corr_t

end