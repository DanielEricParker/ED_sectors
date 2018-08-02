#####################################################################
#This file gives tool for constructing fairly arbitrary spin Hamiltonians

#The main tool is an "abstract Hamiltonian" (i.e. human-readable)
#type called Hamiltonian

#and a way to convert HAMILTONIAN + BASIS into an actual matrix

#### This whole file has gotten super messy and should be cleaned up
#### I think there should be many fewer types

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
    	operator :: OP) = new(prefactor,[operator])

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

    #this is stupidly messy and shows that I need to refactor this code
    ABSTRACT_OP(
    	L :: Int,
    	operator_strng :: String,
    	starting_site :: Int;
    	pbc = false
    	) = ABSTRACT_OP(L,operator_strng,pbc,
    		make_terms_array(L, 
	    		false, 
	    		TERMS(1.0,operator_strng,0,0),
	    		[starting_site])
    	)
end

function Base.show(io::IO, H::ABSTRACT_OP)
	print(io, "ABSTRACT_OP: ",  H.name, " (",length(H.terms), " operators on ", H.L, " sites. PBC = ", H.pbc,")")
	for t in H.terms
		print(io,"\n", t)
	end
end


#ABSTRACT_OP's should form a vector space with
#addition and scalar multiplication

"""
Adds two abstract operators. Operators must be of the same size.
"""
function Base.:+(O1 :: ABSTRACT_OP, O2 :: ABSTRACT_OP)
    if O1.L == O2.L && O1.pbc == O2.pbc
    	newName = string(O1.name," + ",O2.name)
    	newTerms = [O1.terms;O2.terms]
    	return ABSTRACT_OP(O1.L,newName,O1.pbc,newTerms)
    else
    	error("Error ABSTRACT_OP's must have the same length and boundary conditions to be added.")
	end
end

"""
Adds two abstract operators. Operators must be of the same size.
"""
function Base.:*(O :: ABSTRACT_OP, lambda :: Float64)
	newTerms = O.terms
	for t in newTerms
		t.prefactor *= lambda
	end
	return ABSTRACT_OP(O.L,O.name,O.pbc,newTerms)
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


	termsEnd = pbc ? L-1 : L - length(terms.operatorNames)
	sites = collect(terms.startsite : terms.spacing : termsEnd)

	return make_terms_array(L,pbc,terms,sites)
end


function make_terms_array(
	L :: Int,
	pbc :: Bool,
	terms :: TERMS,
	sites :: Array{Int}
	)

	#turn the terms in to abstract terms
	#abstract terms are tuple (prefactor, [Op names])
	allowed_operator_names = ['I','X','Y','Z','+','-','N','C','D','A','B']

	for op_name in terms.operatorNames
		if !in(op_name, allowed_operator_names)
			error("Unsupported operator name \"", op_name, "\" in term, ", terms)
		end
	end 
	opNames = [string(ch) for ch in terms.operatorNames]

	termsArray = []
	for i in sites
		operators = [OP(op_name,(i+j-1)%L) for (j,op_name) in enumerate(opNames)]
		#remove the Identity operators since we don't need those
		operators = filter( (op) -> op.name != "I", operators)
		term = TERM(terms.prefactor, operators)
		push!(termsArray,term)
	end

	return Array{TERM}(termsArray)
end