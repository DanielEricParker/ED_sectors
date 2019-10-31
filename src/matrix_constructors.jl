#####################################################################
#This file gives functions for converting a Basis + ABSTRACT_OP
#into an actual matrix

#####################################################################



"""

    VEC(factor,state)

A type for storing a prefactor * basis_state.
"""
struct VEC
    
    factor :: Complex{Float64}
    state :: UInt64
end
Base.show(io::IO, v::VEC) = print(io, "s[",v.factor, " * ", digits(v.state,base=2), "]")

"""

    is_Null_VEC(vec)

Tests if 'vec' is the zero vector to numerical precision. Mostly used for internal testing.
"""
function is_Null_VEC(v :: VEC)
    #tests if we ended up with zero = VEC(0,0)
    return abs(v.factor) < zeps
end

###################### Operators and operator application ######################

"""

    apply_S_plus(v,i)

Given a basis vector |v>, return S_i^+ |v>, which is another basis vector.

# Arguments
- 'v: VEC': a basis vector VEC(prefactor, state)
- 'i: Int': the index of the spin to apply the operator to, 0 <= i <= L
"""
function apply_S_plus(v::VEC, i :: Int)::VEC
    if ((v.state >> i) & 1 == 0)
        return VEC(2*v.factor, xor(v.state,1 << i))
    else
        return VEC(0,0)
    end 
end


"""

    apply_S_minus(v,i)

#Arguments
*'v: VEC': a basis vector VEC(prefactor, state)
*'i: Int': the index of the spin to apply the operator to, 0 <= i <= L

Given a basis vector |v>, return S_i^- |v>, which is another basis vector.
"""
function apply_S_minus(v::VEC,i::Int)::VEC
    if ((v.state >> i) & 1 == 1)
        return VEC(2*v.factor, xor(v.state,1 << i))
    else
        return VEC(0,0)
    end 
end


"""

    apply_S_z(v,i)

# Arguments
- 'v: VEC': a basis vector VEC(prefactor, state)
- 'i: Int': the index of the spin to apply the operator to, 0 <= i <= L

Given a basis vector |v>, return S_i^z |v>, which is another basis vector.
"""
function apply_S_z(v :: VEC, i :: Int):: VEC
    #x - state
    #i - place to apply the operator 
    if ((v.state >> i) & 1 == 1)
        return VEC(1*v.factor,v.state)
    else
        return VEC(-1*v.factor,v.state)
    end
end

"""

    apply_S_x(v,i)


Given a basis vector |v>, return S_i^x |v>, which is another basis vector.

# Arguments
- 'v: VEC': a basis vector VEC(prefactor, state)
- 'i: Int': the index of the spin to apply the operator to, 0 <= i <= L
"""
function apply_S_x(v :: VEC, i :: Int)::VEC
    return VEC(v.factor,xor(v.state,1 << i))
end


"""

    apply_S_y(v,i)


Given a basis vector |v>, return S_i^y |v>, which is another basis vector.

# Arguments
- 'v: VEC': a basis vector VEC(prefactor, state)
- 'i: Int': the index of the spin to apply the operator to, 0 <= i <= L
"""
function apply_S_y(v :: VEC, i :: Int)::VEC
    flipped = xor(v.state, 1 << i)
    factor = im*v.factor
    if ((v.state >> i) & 1 == 1) #if v[i] == up
        factor *= -1
    end
    return VEC(factor,flipped)
end


"""

    apply_operators(v,t)

Given a basis vector |v> and a term (i.e. operator) t, returns t | v>. For spin-1/2 with the X,Y,Z,+,- operators, the return vector is always another basis vector, rather than superposition thereof, which simplifies things.

# Arguments
- 'v::VEC': a vector VEC(prefactor,basis vector) to apply the term to
- 't::TERM_INTERNAL': a single term from a Hamiltonian to apply

"""
function apply_operators(v :: VEC, t :: TERM_INTERNAL)::VEC
    for op in t.operator
        if is_Null_VEC(v) || op.name == "I" || op.name == "1"
            return v
        else 
            if op.name == "X"
                v = apply_S_x(v,op.site)
            elseif op.name == "Y"
                v = apply_S_y(v,op.site)
            elseif op.name == "Z"
                v = apply_S_z(v,op.site)
            elseif op.name == "+"
                v = apply_S_plus(v,op.site)
            elseif op.name == "-"
                v = apply_S_minus(v,op.site)
            else 
                error("Unsupported operator name: ",op.name)
            end
        end
    end
    v = VEC(v.factor * t.prefactor,v.state)
    return v
end

################ Hamiltonian constructor for a given basis ##########################

"""

    Operator(abstract_op, basis)
    Operator(term, basis)
    Operator(operatorName, site, basis)

Constructs an operator in a given `BASIS`. One can specify a full `ABSTRACT_OP` or, as a shortcut for simple obserables, just a single `TERM` or even the name and site of the operator.

Returns a sparse matrix.



# Examples

Simple observables are simple to create.

```jldoctest
julia> L = 4; basis = BASIS(L);

julia> Sx2 = Operator("X", 2, basis)
16×16 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 16 stored entries:
  [5 ,  1]  =  1.0+0.0im
  [6 ,  2]  =  1.0+0.0im
  [7 ,  3]  =  1.0+0.0im
  [8 ,  4]  =  1.0+0.0im
  [1 ,  5]  =  1.0+0.0im
  [2 ,  6]  =  1.0+0.0im
  [3 ,  7]  =  1.0+0.0im
  [4 ,  8]  =  1.0+0.0im
  [13,  9]  =  1.0+0.0im
  [14, 10]  =  1.0+0.0im
  [15, 11]  =  1.0+0.0im
  [16, 12]  =  1.0+0.0im
  [9 , 13]  =  1.0+0.0im
  [10, 14]  =  1.0+0.0im
  [11, 15]  =  1.0+0.0im
  [12, 16]  =  1.0+0.0im

julia> magnetic_order_parameter = Operator((1/L)*TERM("Z"), basis)
16×16 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 10 stored entries:
  [1 ,  1]  =  1.0+0.0im
  [2 ,  2]  =  0.5+0.0im
  [3 ,  3]  =  0.5+0.0im
  [5 ,  5]  =  0.5+0.0im
  [8 ,  8]  =  -0.5+0.0im
  [9 ,  9]  =  0.5+0.0im
  [12, 12]  =  -0.5+0.0im
  [14, 14]  =  -0.5+0.0im
  [15, 15]  =  -0.5+0.0im
  [16, 16]  =  -1.0+0.0im

```

Creating more involved operators, like most Hamiltonians, is also straightforward.

```jldoctest
julia> L = 4; basis = BASIS(L);

julia> ising = ABSTRACT_OP(L; name="Ising Model", pbc=true) + 2TERM("ZZ") + TERM("X");

julia> H = Operator(ising,basis)
16×16 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 68 stored entries:
  [1 ,  1]  =  8.0+0.0im
  [2 ,  1]  =  1.0+0.0im
  [3 ,  1]  =  1.0+0.0im
  [5 ,  1]  =  1.0+0.0im
  [9 ,  1]  =  1.0+0.0im
  [1 ,  2]  =  1.0+0.0im
  [4 ,  2]  =  1.0+0.0im
  [6 ,  2]  =  1.0+0.0im
  [10,  2]  =  1.0+0.0im
  ⋮
  [7 , 15]  =  1.0+0.0im
  [11, 15]  =  1.0+0.0im
  [13, 15]  =  1.0+0.0im
  [16, 15]  =  1.0+0.0im
  [8 , 16]  =  1.0+0.0im
  [12, 16]  =  1.0+0.0im
  [14, 16]  =  1.0+0.0im
  [15, 16]  =  1.0+0.0im
  [16, 16]  =  8.0+0.0im

```


See also: [`ABSTRACT_OP`](@ref), [`BASIS`](@ref).
"""
function Operator(
    abstract_op :: ABSTRACT_OP,
    basis :: BASIS
    )
    #this should eventually decide if the operator satisfies the symmetries

    #decide if we want the full constructor, or the one with symmetries
    if length(basis.q_numbers) == 0
        construct_matrix_full(basis,abstract_op)
    else 
        construct_matrix_sym(basis,abstract_op)
    end
end

#shortcut constructor for short observables
function Operator(
        term :: TERM,
        basis :: BASIS
    )
    abstract_op  = ABSTRACT_OP(basis.L,term)
    op_mat = Operator(abstract_op, basis)
    return op_mat
end

#shortcut constructor for short observables
function Operator(
        operatorName :: String,
        site :: Int,
        basis :: BASIS
    )
    abstract_op  = ABSTRACT_OP(basis.L,operatorName,site)
    op_mat = Operator(abstract_op, basis)
    return op_mat
end


"""

    construct_matrix_sym(basis, abstract_op)

Internal function to construct a Hamiltonian with symmetries.
"""
function construct_matrix_sym(basis :: BASIS, abstract_op :: ABSTRACT_OP)
    #L - length of the spin chain
    #basis - basis Dict
    #reps - represenatives of conjugacy classes Dict
    #opsList - list of operators we want to add to the Hamiltonian

    dim = basis.dim
    H = spzeros(ComplexF64,dim,dim)

    #iterate over conjugacy classes
    for (x,cc_x) in basis.conj_classes

        #loop over terms in the Hamiltonian
        for term in abstract_op.terms
            #find |y> = H |x> ####ASSUMPTION: this is a state, not superposition
            #this works for X,Y,Z,+,-, but doesn't necessarily hold for everything
            y = apply_operators(VEC(1.0,x),term)

            #check if the element exists --- it could have zero norm
            if ~is_Null_VEC(y) && haskey(basis.get_conj_class,y.state)

                bv_y = basis.get_conj_class[y.state]
                cc_y = basis.conj_classes[bv_y.conj_class]

                #check if upper diagonal
                if cc_y.index >= cc_x.index

                    #println("x: "*string(digits(x,2,4))*", cc: $(cc_x), y: "*string(digits(y.state,2,4))*", cc: $(cc_y), factor: $(y.factor)")
                    #compute the matrix element --- see Note
                    h_xy = zchop(sqrt(cc_y.norm/cc_x.norm) * bv_y.phase_factor * y.factor)
                    H[cc_x.index,cc_y.index] += h_xy

                    #check if diagonal
                    if cc_y.index > cc_x.index
                         H[cc_y.index,cc_x.index] += conj(h_xy)
                    end
                end
                # else
                #     println("Warning: no target in basis for state ",y.state)
                #end
            end
        end
    end
    return H
    #Making it Hermitian makes diagonalization ridiculously slow for unclear reasons right now
    # #maybe a bug?
    ## see https://discourse.julialang.org/t/unexpectedly-slow-performance-of-eigs-with-a-hermitian-view/13701

    # return Hermitian(H,:U)
end 

###define Pauli operators concretely

const pauli_1 = sparse([1,2],[1,2],[1.0,1.0])
const pauli_x = sparse([1,2],[2,1],[1.0,1.0])
const pauli_y = sparse([1,2],[2,1],[1.0*im,-1.0*im])
const pauli_z = sparse([1,2],[1,2],[1.0,-1.0])


const pauli_plus = sparse([0 2.0; 0 0])
const pauli_minus = sparse([0 0; 2.0 0])


# ###also do spinless fermions

# const cdag_op = sparse([0 2.0; 0 0])
# const c_op = sparse([0 0; 2.0 0])
# const number_op = sparse([1.0 0.0; 0.0 0.0])


 
"""
Returns the correct Pauli matrix for the name. Really this should a hardcoded dictionary, but its essentially the same and there's not that many of them.
"""
function get_pauli_matrix(name :: String)
        if name == "X"
            return pauli_x
        elseif name == "Y"
            return pauli_y
        elseif name == "Z"
            return pauli_z
        elseif name == "+"
            return pauli_plus
        elseif name == "-"
            return pauli_minus
        elseif name == "I" || name == "1"
            return pauli_1
        # elseif name == "C" #destroy spinless fermions
        #     return c_op
        # elseif name == "D" #C dagger
        #     return cdag_op
        # elseif name == "N"
        #     return number_op
        else 
            error("Unsupported operator name: ",name)
        end
end

#formally untested, but this has been verified against other code from matlab, and a few other things
"""

    construct_matrix_full(basis, abstract_op)

An internal method to quickly make an operator for the full basis with no symmetry constraints. 
"""
function construct_matrix_full(basis :: BASIS, abstract_op :: ABSTRACT_OP)

    dim = 2^basis.L
    H = spzeros(ComplexF64,dim,dim)

    #loop over terms in the Hamiltonian
    for term in abstract_op.terms

        #start with the Identity
        term_matrix = sparse(I,1,1)
        last_pos = -1

        #loop over terms, i.e. 2.0 XZX -> [X,Z,X]
        for op in term.operator
            pos = op.site
            #how many ones in the middle?
            pos_shift = pos-last_pos-1
            #get the pauli matrix
            site_mat = get_pauli_matrix(op.name)

            #tensor product
            if pos_shift == 0
                term_matrix = kron(site_mat,term_matrix)
            else 
                shift_size = 2^pos_shift
                term_matrix = kron(kron(site_mat,sparse(I, shift_size, shift_size)),term_matrix)
            end
            last_pos = pos
        end

        #fill out matrix to 2^L x 2^L
        pos_shift = basis.L - 1 - last_pos
        if pos_shift > 0
            shift_size = 2^pos_shift
            term_matrix = kron(sparse(I,shift_size,shift_size),term_matrix)
        end
        
        #multiply by the prefactor
        term_matrix *= term.prefactor

        #add the term
        H += term_matrix
    end

    return H
end

#######################################
# New and more efficient version of the constructor
#######################################



"""
Efficiently applies a pauli string on a state
"""
function Base.:*(P :: PAULI, psi :: VEC) :: VEC
	pre = P.x * psi.factor

	state_new = xor(P.V.w,psi.state)
	minus_ones = count_ones(P.V.v & (~state_new)) % 2 == 1

	if xor(P.V.epsilon, minus_ones)
		pre = -pre
	end

	if P.V.delta
		pre = Complex(-pre.im,pre.re)
	end


	return VEC(pre, state_new)
end


"""
Efficiently applies a pauli string on a state
"""
function Base.:*(P :: PAULI, x :: UInt64) :: VEC
	pre = P.x
	state_new = xor(P.V.w,x)
	minus_ones = count_ones(P.V.v & (~state_new)) % 2 == 1

	if xor(P.V.epsilon, minus_ones)
		pre = -pre
	end
	if P.V.delta
		pre = Complex(-pre.im,pre.re)
	end

	return VEC(pre, state_new)
end




"""
Converts a TERM_INTERNAL to a PAULI string.
"""
#This is currently done in a terribly stupid way.
function PauliString(term :: TERM_INTERNAL) :: PAULI

	max_site = maximum([o.site for o in term.operator])
	# println(max_site)

	str_array = ["I" for i in 0:max_site]
    # println(str_array)
	for o in term.operator
		str_array[o.site+1] = o.name
	end
	pauli_string = join(reverse(str_array))
	# println("string: ",pauli_string)

	return term.prefactor * PAULI(pauli_string)
end

"""

"""




#formally untested, but this has been verified against other code from matlab, and a few other things
"""

    construct_matrix_full(basis, abstract_op)

An internal method to quickly make an operator for the full basis with no symmetry constraints. 
"""
function construct_matrix_full_new(basis :: BASIS, abstract_op :: ABSTRACT_OP)

    P_op = PauliOp(abstract_op)

    return construct_matrix_full_new_Pauli(basis,P_op)
end

function construct_matrix_full_new_Pauli(basis :: BASIS, PauliStrings ::PauliOp)

    dim = 2^basis.L
    H = spzeros(ComplexF64,dim,dim)

    #loop over terms in the Hamiltonian
    for n in 0:dim-1
      for (V,x) in PauliStrings
        m = PAULI(x,V)*UInt64(n)
        # println(P,"*",x,"=",y)
        H[dim-m.state,dim-n] += m.factor
      end
  end

    return H
end


function PauliOp(abstract_op :: ABSTRACT_OP) :: PauliOp
  O = Dict{VERT,Complex{Float64}}()
  for term in consolidate_coefficients(abstract_op).terms
    P = PauliString(term)

    if haskey(O, P.V)
      O[P.V] += P.x
    else
      O[P.V] = P.x
    end
  end

  return O
end



# #UNTESTED!
# """

#     construct_matrix_fermion_biliners(basis, abstract_op)

# Makes a Hamiltonian from fermion bilinears. Only work for number-conserving Hamiltonians

# # Argument s
# - 'basis :: Basis': the basis for the spin chain
# - 'abstract_op :: HAMILTONIAN': the abstract operator to implement
# """
# function construct_matrix_fermion_biliners(basis :: BASIS, abstract_op :: ABSTRACT_OP)
#     dim = basis.L
#     H = spzeros(ComplexF64,dim,dim)

#     #loop over terms in the Hamiltonian
#     for term in abstract_op.terms
#         op = term.operator
#         #is it the number operator?
#         if length(op) == 1 && op[1].name == "N"
#             i = op[1].site
#             j = i
#         else
#             #check ordering
#             if op[1].name == "D" && op[2].name == "C"
#                 i = op[1].site
#                 j = op[2].site
#             elseif op[1].name == "C" && op[2].name == "D"
#                 i = op[2].site
#                 j = op[1].site
#             else
#                 error("Unknown term:", term)
#             end
#         end
#         # println(term)
#         # println("H[",i,",",j, "] = ", term.prefactor)
#         H[i+1,j+1]  = term.prefactor
#     end

#     return H
# end

#UNTESTED!

# """
# Makes a Hamiltonian for a single mode of bosons.
# #Argument 
# * 'N :: Int': the maximum number of modes possible
# * 'abstract_op :: HAMILTONIAN': the abstract operator to implement
# """

# function construct_matrix_boson_mode(basis :: Basis, abstract_op :: ABSTRACT_OP)

#     dim = N
#     H = spzeros(ComplexF64, dim, dim)

#     #loop over terms in the Hamiltonian
#     for term in abstract_op.terms
#         op = term.operator

#         if length(op) == 1 && op[1].name == "I"
#             #this is silly, but I don't understand how sparse matrix constructors work
#             H += term.prefactor*sparse(1:N,1:N,1.0)
#         elseif length(op) == 1 && op[1].name == "A"
#             #raising operator
#             H += term.prefactor*sparse(1:N-1,2:N,1.0,N,N)
#         elseif length(op) == 1 && op[1].name == "B"
#             #lowering operator
#             H += term.prefactor*sparse(1:N-1,2:N,1.0,N,N)
#         elseif length(op) == 1 && op[1].name == "N"
#             #number operator
#             H += term.prefactor*sparse(1:N,1:N,N:-1:1,N,N)
#         else
#             error("Unknown term:", term)
#         end
#     return H
# end