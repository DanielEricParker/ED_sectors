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

	EIGSYS(evs, U)
	EIGSYS(H)

A struct for storing the full spectrum `evs` and eigenstates of a system `U`. Explicitly, `U` is the (unitary) matrix which gives the change of basis from the original to energy bases, i.e. the *columns* of `U` are the eigenvectors of the system. So for an `EIGSYS` `ES`, the full spectrum is `ES.evs` and the 3rd eigenvector, for example, is `ES.U[:,3]`.

# Examples
```jldoctest
julia> L = 4; basis = BASIS(L);

julia> ising = ABSTRACT_OP(L; name="Ising Model", pbc=true) + 2TERM("ZZ") + TERM("X");

julia> H = Operator(ising,basis);

julia> ES = EIGSYS(H)
EIGSYS[dim: 16]

julia> psi3 = ES.U[:,3]
16-element Array{Complex{Float64},1}:
  0.08727992975105786 + 0.0im
 -0.23235256582864727 + 0.0im
  -0.2323525658286471 + 0.0im
   0.3509044118546869 + 0.0im
 -0.23235256582864688 + 0.0im
 -0.17367654405779392 + 0.0im
  0.35090441185468646 + 0.0im
 -0.23235256582864808 + 0.0im
 -0.23235256582864716 + 0.0im
   0.3509044118546871 + 0.0im
  -0.1736765440577938 + 0.0im
 -0.23235256582864866 + 0.0im
  0.35090441185468646 + 0.0im
 -0.23235256582864844 + 0.0im
 -0.23235256582864783 + 0.0im
  0.08727992975105792 - 0.0im

julia> isapprox((H * psi3),(ES.evs[3]*psi3)) #they are the same to numerical error
true

```
"""
struct EIGSYS
	evs :: Array{Float64} #eigenvalues
	U :: Array{Complex{Float64},2} #eigenvectors


	"""
	Function to perform full diagonalization on the Hamiltonian.
	Returns an EigSys.
	"""
	function EIGSYS(H :: AbstractMatrix)
		eigensystem = eigen(Matrix(H))
		new(real(eigensystem.values),eigensystem.vectors)
	end


	"""
	Alternatively, construct some other way.
	Returns an EigSys.
	"""
	function EIGSYS(
				evs :: Array{Float64}, #eigenvalues
				U :: Array{Complex{Float64},2} #eigenvectors
			)
		new(evs,O)
	end
end
Base.show(io::IO, eigsys :: EIGSYS) = print(io, "EIGSYS[dim: $(length(eigsys.evs))]")


# """
# Function to compute n-point correlation functions for a density matrix rho.

# """
# function correlation_npt(
# 	rho :: Array{Complex{Float64},2},
# 	Ops :: Array{SOP},
# 	factor :: Float64,
# 	basis :: BASIS
# 	)

# end

"""

	thermal_density_matrix(eigsys, beta)

Computes the thermal density matrix *in the energy basis* for a given `eigsys` and inverse temperature `beta`. When `beta = 0`, this gives the "infinite-temperature" density matrix.

!!! warning

    It is easy to get underflow or overflow errors by exponentiating large numbers. The largest Float64 is roughly ``10^{308} \\sim e^{709}``. It is recommended to double-check that your results make sense after using this function.

"""
function thermal_density_matrix(
	eigsys :: EIGSYS,
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
	eigsys :: EIGSYS,
	beta :: Int
	)
	if(beta != 0)
		return thermal_density_matrix(eigsys,Float64(beta))
	end
	D = Matrix{ComplexF64}(Diagonal([1/length(eigsys.evs) for i in 1:Z]))
	return D
end

#verified against mathematica code
#this is easy to read rather than fast, but it shouldn't matter much
"""

	r_statistic(eigenvalues :: Array{Float64})
	r_statistic(ES :: EIGSYS)

Computes the R-statistic for a spectrum. The R-statistic is a measure of how close to integrability the system is. For integrable systems, ``\\left<r\\right> \\sim 0.386``. For GOE level-statistics, ``\\left<r\\right> \\sim 0.5295``. 

!!! warning

    Using the R statistic correctly is subtle. It only makes sense systems sufficiently large as to be close to the thermodynamic limit. Moreover, one must restrict to a single symmetry sector. Any residual symmetry will give an artificially low R-statistic. See, e.g., <https://arxiv.org/abs/cond-mat/0610854> for details..

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


function r_statistic(
	eigsys :: EIGSYS
	)
	return r_statistic(eigsys.evs)
end