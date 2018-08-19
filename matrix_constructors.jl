#####################################################################
#This file gives functions for converting a Basis + ABSTRACT_OP
#into an actual matrix

#####################################################################



###################### Operators and operator application ######################

"""

    apply_S_plus(v,i)

#Arguments
*'v: VEC': a basis vector VEC(prefactor, state)
*'i: Int': the index of the spin to apply the operator to, 0 <= i <= L

Given a basis vector |v>, return S_i^+ |v>, which is another basis vector.
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

#Arguments
*'v: VEC': a basis vector VEC(prefactor, state)
*'i: Int': the index of the spin to apply the operator to, 0 <= i <= L

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

#Arguments
*'v: VEC': a basis vector VEC(prefactor, state)
*'i: Int': the index of the spin to apply the operator to, 0 <= i <= L

Given a basis vector |v>, return S_i^x |v>, which is another basis vector.
"""
function apply_S_x(v :: VEC, i :: Int)::VEC
    return VEC(v.factor,xor(v.state,1 << i))
end


"""

    apply_S_y(v,i)

#Arguments
*'v: VEC': a basis vector VEC(prefactor, state)
*'i: Int': the index of the spin to apply the operator to, 0 <= i <= L

Given a basis vector |v>, return S_i^y |v>, which is another basis vector.
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

#Arguments
*'v::VEC': a vector VEC(prefactor,basis vector) to apply the term to
*'t::TERM': a single term from a Hamiltonian to apply

Given a basis vector |v> and a term (i.e. operator) t, returns t | v>. For spin-1/2 with the X,Y,Z,+,- operators, the return vector is always another basis vector, rather than superposition thereof, which simplifies things.
"""
function apply_operators(v :: VEC, t :: TERM)::VEC
    for op in t.operator
        if is_Null_VEC(v)
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

    construct_matrix(basis, abstract_op)

Constructs a matrix from an abstract operator for the given basis. Returns a sparse matrix.
"""
function construct_matrix(
    basis :: Basis, 
    abstract_op :: ABSTRACT_OP)

    #this should eventually decide if the operator satisfies the symmetries

    #decide if we want the full constructor, or the one with symmetries
    if length(basis.q_numbers) == 0
        construct_matrix_full(basis,abstract_op)
    else 
        construct_matrix_sym(basis,abstract_op)
    end
end

"""

    construct_matrix_sym(basis, abstract_op)

Internal function to construct a Hamiltonian with symmetries.
"""
function construct_matrix_sym(basis :: Basis, abstract_op :: ABSTRACT_OP)
    #L - length of the spin chain
    #basis - basis Dict
    #reps - represenatives of conjugacy classes Dict
    #opsList - list of operators we want to add to the Hamiltonian

    dim = length(basis.conj_classes)
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

if testing2
    println("Testing making the Hamiltonian from the full basis")
    abstract_hamiltonian = make_XXZ_operators_new(4,0.8)
   	#println(abstract_hamiltonian)
    H2 = make_Hamiltonian(4,basisFull,abstract_hamiltonian)
    println(H1 == H2)
end




###define Pauli operators concretely

const pauli_1 = sparse([1,2],[1,2],[1.0,1.0])
const pauli_x = sparse([1,2],[2,1],[1.0,1.0])
const pauli_y = sparse([1,2],[2,1],[1.0*im,-1.0*im])
const pauli_z = sparse([1,2],[1,2],[1.0,-1.0])


const pauli_plus = sparse([0 2.0; 0 0])
const pauli_minus = sparse([0 0; 2.0 0])


###also do spinless fermions

const cdag_op = sparse([0 2.0; 0 0])
const c_op = sparse([0 0; 2.0 0])
const number_op = sparse([1.0 0.0; 0.0 0.0])


 
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
        elseif name == "C" #destroy spinless fermions
            return c_op
        elseif name == "D" #C dagger
            return cdag_op
        elseif name == "N"
            return number_op
        else 
            error("Unsupported operator name: ",name)
        end
end


"""

    construct_matrix_full(basis, abstract_op)

An internal method to quickly make an operator for the full basis with no symmetry constraints. 
"""
function construct_matrix_full(basis :: Basis, abstract_op :: ABSTRACT_OP)

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



#UNTESTED!
"""

    construct_matrix_fermion_biliners(basis, abstract_op)

Makes a Hamiltonian from fermion bilinears. Only work for number-conserving Hamiltonians

#Argument 
* 'basis :: Baiss': the basis for the spin chain
* 'abstract_op :: HAMILTONIAN': the abstract operator to implement
"""
function construct_matrix_fermion_biliners(basis :: Basis, abstract_op :: ABSTRACT_OP)
    dim = basis.L
    H = spzeros(ComplexF64,dim,dim)

    #loop over terms in the Hamiltonian
    for term in abstract_op.terms
        op = term.operator
        #is it the number operator?
        if length(op) == 1 && op[1].name == "N"
            i = op[1].site
            j = i
        else
            #check ordering
            if op[1].name == "D" && op[2].name == "C"
                i = op[1].site
                j = op[2].site
            elseif op[1].name == "C" && op[2].name == "D"
                i = op[2].site
                j = op[1].site
            else
                error("Unknown term:", term)
            end
        end
        # println(term)
        # println("H[",i,",",j, "] = ", term.prefactor)
        H[i+1,j+1]  = term.prefactor
    end

    return H
end

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