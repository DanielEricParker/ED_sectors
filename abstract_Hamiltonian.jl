#####################################################################
#This file gives tool for constructing fairly arbitrary spin Hamiltonians

#The main tool is an "abstract Hamiltonian" (i.e. human-readable)
#type called Hamiltonian

#and a way to convert HAMILTONIAN + BASIS into an actual matrix

#####################################################################



#################### Data Structures ############################


struct VEC
    #type for storing factor * basis_state
    factor :: Complex{Float64}
    state :: UInt64
end
Base.show(io::IO, v::VEC) = print(io, "s[",v.factor, " * ", bin(v.state), "]")

if testing
    #TEST VEC
    println(VEC(1.3,7))
end

function is_Null_VEC(v :: VEC)
    #tests if we ended up with zero = VEC(0,0)
    return abs(v.factor) < 1e-15 
end




struct OP
    #type for storing factor * operator_name
    #such as "S_z", "S_x", "S_+"
    name :: String
    site :: Int
end
Base.show(io::IO, o::OP) = print(io,  o.name, "_", o.site)

if testing
    #TEST OP
    println(OP("X",3))
end







struct TERM
    #type for storing a term in a Hamiltonian
    #such as "0.5 Sz Sz"
    prefactor :: Float64
    operator :: Array{OP}

    TERM(prefactor :: Float64,
         operator :: Array{OP}) = new(
             prefactor,
             #we want to automatically sort the operators by site
             sort(operator, by = (op) -> op.site)
         )
end
Base.show(io::IO, t::TERM) = print(io,  t.prefactor, "*", t.operator)

if testing2
    testTerm = TERM(0.5,[OP("+",3),OP("-",4),OP("-",5)])
    println(testTerm)
    println(typeof(testTerm.operator))
end


struct TERMS
	#type for storing a periodically repeated set of terms
	prefactor :: Float64
	operatorNames :: String
	startsite :: Int
	spacing :: Int

	#constructors
	TERMS(prefactor :: Real,
		operatorNames :: String,
		startsite :: Int,
		spacing :: Int) = new(Float64(prefactor), operatorNames, startsite, spacing)

	TERMS(prefactor :: Real, 
		operatorNames :: String) = new(Float64(prefactor), operatorNames, 0, 1)
end
Base.show(io::IO, t::TERMS) = print(io,  t.prefactor, "*", t.operatorNames, " starting at site ", t.startsite, ", spaced by ", t.spacing)


if testing2
	testTerms = TERMS(3.0,"ZXZ",1,2)
	println(testTerms)


	testTerms = TERMS(-4,"ZWXZ")
	println(testTerms)
end




struct ABSTRACT_OP
    #type for storing Hamiltonians
    L :: Int
    name :: String
    pbc :: Bool
    terms :: Array{TERM}



    function ABSTRACT_OP(    
    	L :: Int,
	    name :: String,
	    pbc :: Bool,
	    terms :: Array{TERM})

	    #implements check_Hamiltonian
	    for term in terms
	        for op in term.operator
	            if op.site < 0 || op.site >= L
	                error("Out of bounds! Term:", term)
	            end
	        end
	    end

	    #removed terms that are basically zero for float errors
	    #this should be improved at some point
	    terms = filter( t -> abs(t.prefactor) >= 10e-12, terms)

	    new(L,name,pbc,terms)
	end

    #short constructor for making this quickly
    ABSTRACT_OP(
        L :: Int,
        op :: OP,
        ) = new(L,"",false,[TERM(1.0,[op])])

    ABSTRACT_OP(    
    	L :: Int,
	    name :: String,
	    pbc :: Bool) = new(L,name,pbc,[])


    #non-periodic boundary conditions by default
    ABSTRACT_OP(
    	L :: Int,
    	name :: String,
    	terms :: Array{TERM})= ABSTRACT_OP(L,name,false,terms)
end

function Base.show(io::IO, H::ABSTRACT_OP)
	print(io, "ABSTRACT_OP: ",  H.name, " (",length(H.terms), " operators on ", H.L, " sites. PBC = ", H.pbc,")")
	for t in H.terms
		print(io,"\n", t)
	end
end


"""
Gives a way to add single terms to Hamiltonians
usage: H += TERM(0.5,[OP("+",3),OP("-",4),OP("-",5)])
"""
function Base.:+(H :: ABSTRACT_OP, t :: TERM)
	newTerms = [H.terms; t]    
    return ABSTRACT_OP(H.L,H.name,H.pbc,newTerms)
end

if testing2
    println("Testing making a Hamiltonian")
    t1 = term(0.5,OP("+",3),OP("-",4),OP("-",5))
    t2 = term(0.5,OP("-",3),OP("z",5))
    println(typeof(t1.operator))
    println(typeof(t2.operator))
    testHam = ABSTRACT_OP("testing", [t1])
    println(testHam)
    
    testHam += t2
    println(testHam)
end

"""
Gives a way to add many repeated terms to Hamiltonians
usage: H += TERMS(0.5,"ZZ")
or     H += TERMS(0.5,"ZZ",1,2) #start at site 1, space by 2
"""
function Base.:+(H :: ABSTRACT_OP, terms :: TERMS)

	newTerms = [H.terms; make_terms(H.L,H.pbc,terms)]
	return ABSTRACT_OP(H.L,H.name,H.pbc,newTerms)
end


"""
Parses a TERMS struct into individual TERM structs
"""
function make_terms(
	L :: Int,
	pbc :: Bool,
	terms :: TERMS)

	#turn the terms in to abstract terms
	#abstract terms are tuple (prefactor, [Op names])
	allowed_operator_names = ['I','X','Y','Z','+','-']

	for op_name in terms.operatorNames
		if !in(op_name, allowed_operator_names)
			error("Unsupported operator name \"", op_name, "\" in term, ", terms)
		end
	end
	opNames = [string(ch) for ch in terms.operatorNames]

	termEnd = pbc ? L-1 : L - length(opNames)

	termsArray = []
	for i in terms.startsite : terms.spacing : termEnd
		operators = [OP(op_name,(i+j-1)%L) for (j,op_name) in enumerate(opNames)]
		#remove the Identity operators since we don't need those
		operators = filter( (op) -> op.name != "I", operators)
		term = TERM(terms.prefactor, operators)
		push!(termsArray,term)
	end

	return Array{TERM}(termsArray)
end



###################### Operators and operator application ######################


function apply_S_plus(v::VEC, i :: Int)
    #x - state
    #i - place to apply the operator 
    if ((v.state >> i) & 1 == 0)
        return VEC(2*v.factor, xor(v.state,1 << i))
    else
        return VEC(0,0)
    end 
end

if testing

    #TEST apply_S_plus
    v = VEC(3.1,5)
    w = VEC(3.1,6)
    println("Testing S plus")
    println(v)
    println(apply_S_plus(v,1))
    println(w)
    println(apply_S_plus(w,1))
    println(is_Null_VEC(apply_S_plus(w,1)))
end

function apply_S_minus(v::VEC,i::Int)
    #x - state
    #i - place to apply the operator 
    if ((v.state >> i) & 1 == 1)
        return VEC(2*v.factor, xor(v.state,1 << i))
    else
        return VEC(0,0)
    end 
end

if testing
    println("Testing S minus")
    println(v)
    println(apply_S_minus(v,1))
    println(w)
    println(apply_S_minus(w,1))
    println(is_Null_VEC(apply_S_minus(w,1)))
end


function apply_S_z(v :: VEC, i :: Int)
    #x - state
    #i - place to apply the operator 
    if ((v.state >> i) & 1 == 1)
        return VEC(1*v.factor,v.state)
    else
        return VEC(-1*v.factor,v.state)
    end
end

if testing
    println("Testing S_z")
    println(v)
    println(apply_S_z(v,1))
    println(w)
    println(apply_S_z(w,1))
    println(is_Null_VEC(apply_S_z(w,1)))
end



function apply_S_x(v :: VEC, i :: Int)
    #x - state
    #i - place to apply the operator 
    return VEC(v.factor,xor(v.state,1 << i))
end

if testing
    println("Testing S_x")
    println(v)
    println(apply_S_x(v,1))
    println(w)
    println(apply_S_x(w,1))
    println(is_Null_VEC(apply_S_x(w,1)))
end


####not verified
function apply_S_y(v :: VEC, i :: Int)
    #v - Vector
    #i - place to apply the operator 
    flipped = xor(v.state, 1 << i)
    factor = im*v.factor
    if ((v.state >> i) & 1 == 1) #if v[i] == up
        factor *= -1
    end
    return VEC(factor,flipped)
end

if testing
    println("Testing S_y")
    println(v)
    println(apply_S_y(v,1))
    println(w)
    println(apply_S_y(w,1))
    println(is_Null_VEC(apply_S_y(w,1)))
end






function apply_operators(v :: VEC, t :: TERM)
    # v vector
    # list of operators to apply to various sites
    # only works for operators whose output is a single basis vector
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

if testing2
    ops = [OP("X",1),OP("X",2)]
    ops4 = [OP("-",5)]
    v = VEC(1.0,15)
    println("Testing Applying Multiple Operators")
    println(v)
    map(println,ops)
    println(apply_operators(v,TERM(2.0,ops)))
    println(apply_operators(v,TERM(3.0,ops4)))
end

################ Hamiltonian constructor for a given basis ##########################

"""
Construts a matrix from an abstract operator for the given basis.
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
            y = apply_operators(VEC(1.0,x),term)

            #check if the element exists --- it could have zero norm
            if ~is_Null_VEC(y) && haskey(basis.get_conj_class,y.state)

                bv_y = basis.get_conj_class[y.state]
                cc_y = basis.conj_classes[bv_y.conj_class]

                #check if upper diagonal
                if cc_y.index >= cc_x.index
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

pauli_1 = sparse([1,2],[1,2],[1.0,1.0])
pauli_x = sparse([1,2],[2,1],[1.0,1.0])
pauli_y = sparse([1,2],[2,1],[1.0*im,-1.0*im])
pauli_z = sparse([1,2],[1,2],[1.0,-1.0])


pauli_plus = sparse([0 2.0; 0 0])
pauli_minus = sparse([0 0; 2.0 0])

 
"""
Returns the correct Pauli matrix for the name
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
        else 
            error("Unsupported operator name: ",name)
        end
end




"""
Quickly makes a Hamiltonian for the full basis with no symmetry constraints. 

#Argument 
* 'L :: Int': the length of the spin chain
* 'abstract_Ham :: HAMILTONIAN': the abstract operator to implement
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
                term_matrix = kron(term_matrix,site_mat)
            else 
                shift_size = 2^pos_shift
                term_matrix = kron(term_matrix,kron(sparse(I, shift_size, shift_size),site_mat))
            end
            last_pos = pos
        end

        #fill out matrix to 2^L x 2^L
        pos_shift = basis.L - 1 - last_pos
        if pos_shift > 0
            shift_size = 2^pos_shift
            term_matrix = kron(term_matrix, sparse(I,shift_size,shift_size))
        end
        
        #multiply by the prefactor
        term_matrix *= term.prefactor

        #add the term
        H += term_matrix
    end

    return H
end






