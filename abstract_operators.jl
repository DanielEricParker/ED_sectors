#####################################################################
#This file gives tool for constructing fairly arbitrary spin Hamiltonians

#The main tool is an "abstract Hamiltonian" (i.e. human-readable)
#type called Hamiltonian

#and a way to convert HAMILTONIAN + BASIS into an actual matrix

#### This whole file has gotten super messy and should be cleaned up
#### I think there should be many fewer types

#####################################################################




# #################### Data Structures ############################

# struct OP
#     #type for storing factor * operator_name
#     #such as "S_z", "S_x", "S_+"
#     name :: String
#     site :: Int
# end
# Base.show(io::IO, o::OP) = print(io,  o.name, "_", o.site)

# if testing
#     #TEST OP
#     println(OP("X",3))
# end

# struct TERM
#     #type for storing a term in a Hamiltonian
#     #such as "0.5 Sz Sz"
#     prefactor :: Float64
#     operator :: Array{OP}

#     TERM(prefactor :: Float64,
#     	operator :: OP) = new(prefactor,[operator])

#     TERM(prefactor :: Float64,
#          operator :: Array{OP}) = new(
#              prefactor,
#              #we want to automatically sort the operators by site
#              sort(operator, by = (op) -> op.site)
#          )
# end
# Base.show(io::IO, t::TERM) = print(io,  t.prefactor, "*", t.operator)

# if testing2
#     testTerm = TERM(0.5,[OP("+",3),OP("-",4),OP("-",5)])
#     println(testTerm)
# 	println(typeof(testTerm.operator))
# end


# struct TERMS
# 	#type for storing a periodically repeated set of terms
# 	prefactor :: Float64
# 	operatorNames :: String
# 	startsite :: Int
# 	spacing :: Int

# 	#constructors
# 	TERMS(prefactor :: Real,
# 		operatorNames :: String,
# 		startsite :: Int,
# 		spacing :: Int) = new(Float64(prefactor), operatorNames, startsite, spacing)

# 	TERMS(prefactor :: Real, 
# 		operatorNames :: String) = new(Float64(prefactor), operatorNames, 0, 1)
# end
# Base.show(io::IO, t::TERMS) = print(io,  t.prefactor, "*", t.operatorNames, " starting at site ", t.startsite, ", spaced by ", t.spacing)


# if testing2
# 	testTerms = TERMS(3.0,"ZXZ",1,2)
# 	println(testTerms)


# 	testTerms = TERMS(-4,"ZWXZ")
# 	println(testTerms)
# end


# struct ABSTRACT_OP
#     #type for storing Hamiltonians
#     L :: Int
#     name :: String
#     pbc :: Bool
#     terms :: Array{TERM}



#     function ABSTRACT_OP(    
#     	L :: Int,
# 	    name :: String,
# 	    pbc :: Bool,
# 	    terms :: Array{TERM})

# 	    #implements check_Hamiltonian
# 	    for term in terms
# 	        for op in term.operator
# 	            if op.site < 0 || op.site >= L
# 	                error("Out of bounds! Term:", term)
# 	            end
# 	        end
# 	    end

# 	    #removed terms that are basically zero for float errors
# 	    #this should be improved at some point
# 	    terms = filter( t -> abs(t.prefactor) >= 10e-12, terms)

# 	    new(L,name,pbc,terms)
# 	end

#     #short constructor for making this quickly
#     ABSTRACT_OP(
#         L :: Int,
#         op :: OP,
#         ) = new(L,"",false,[TERM(1.0,[op])])

#     ABSTRACT_OP(    
#     	L :: Int,
# 	    name :: String,
# 	    pbc :: Bool) = new(L,name,pbc,[])


#     #non-periodic boundary conditions by default
#     ABSTRACT_OP(
#     	L :: Int,
#     	name :: String,
#     	terms :: Array{TERM})= ABSTRACT_OP(L,name,false,terms)

#     #this is stupidly messy and shows that I need to refactor this code
#     ABSTRACT_OP(
#     	L :: Int, 
#     	operator_strng :: String,
#     	starting_site :: Int;
#     	pbc = false
#     	) = ABSTRACT_OP(L,operator_strng,pbc,
#     		make_terms_array(L, 
# 	    		false, 
# 	    		TERMS(1.0,operator_strng,0,0),
# 	    		[starting_site])
#     	)
# end

# function Base.show(io::IO, H::ABSTRACT_OP)
# 	print(io, "ABSTRACT_OP: ",  H.name, " (",length(H.terms), " operators on ", H.L, " sites. PBC = ", H.pbc,")")
# 	for t in H.terms
# 		print(io,"\n", t)
# 	end
# end


# #ABSTRACT_OP's should form a vector space with
# #addition and scalar multiplication

# """
# Adds two abstract operators. Operators must be of the same size.
# """
# function Base.:+(O1 :: ABSTRACT_OP, O2 :: ABSTRACT_OP)
#     if O1.L == O2.L && O1.pbc == O2.pbc
#     	newName = string(O1.name," + ",O2.name)
#     	newTerms = [O1.terms;O2.terms]
#     	return ABSTRACT_OP(O1.L,newName,O1.pbc,newTerms)
#     else
#     	error("Error ABSTRACT_OP's must have the same length and boundary conditions to be added.")
# 	end
# end

# """
# Adds two abstract operators. Operators must be of the same size.
# """
# function Base.:*(O :: ABSTRACT_OP, lambda :: Float64)
# 	newTerms = O.terms
# 	for t in newTerms
# 		t.prefactor *= lambda
# 	end
# 	return ABSTRACT_OP(O.L,O.name,O.pbc,newTerms)
# end



# """
# Gives a way to add single terms to Hamiltonians
# usage: H += TERM(0.5,[OP("+",3),OP("-",4),OP("-",5)])
# """
# function Base.:+(H :: ABSTRACT_OP, t :: TERM)
# 	newTerms = [H.terms; t]    
#     return ABSTRACT_OP(H.L,H.name,H.pbc,newTerms)
# end

# if testing2
#     println("Testing making a Hamiltonian")
#     t1 = term(0.5,OP("+",3),OP("-",4),OP("-",5))
#     t2 = term(0.5,OP("-",3),OP("z",5))
#     println(typeof(t1.operator))
#     println(typeof(t2.operator))
#     testHam = ABSTRACT_OP("testing", [t1])
#     println(testHam)
    
#     testHam += t2
#     println(testHam)
# end

# """
# Gives a way to add many repeated terms to Hamiltonians
# usage: H += TERMS(0.5,"ZZ")
# or     H += TERMS(0.5,"ZZ",1,2) #start at site 1, space by 2
# """
# function Base.:+(H :: ABSTRACT_OP, terms :: TERMS)

# 	newTerms = [H.terms; make_terms(H.L,H.pbc,terms)]
# 	return ABSTRACT_OP(H.L,H.name,H.pbc,newTerms)
# end



# """
# Parses a TERMS struct into individual TERM structs
# """
# function make_terms(
# 	L :: Int,
# 	pbc :: Bool,
# 	terms :: TERMS)


# 	termsEnd = pbc ? L-1 : L - length(terms.operatorNames)
# 	sites = collect(terms.startsite : terms.spacing : termsEnd)

# 	return make_terms_array(L,pbc,terms,sites)
# end


# function make_terms_array(
# 	L :: Int,
# 	pbc :: Bool,
# 	terms :: TERMS,
# 	sites :: Array{Int}
# 	)

# 	#turn the terms in to abstract terms
# 	#abstract terms are tuple (prefactor, [Op names])
# 	allowed_operator_names = ['I','X','Y','Z','+','-','N','C','D','A','B']

# 	for op_name in terms.operatorNames
# 		if !in(op_name, allowed_operator_names)
# 			error("Unsupported operator name \"", op_name, "\" in term, ", terms)
# 		end
# 	end 
# 	opNames = [string(ch) for ch in terms.operatorNames]

# 	termsArray = []
# 	for i in sites
# 		operators = [OP(op_name,(i+j-1)%L) for (j,op_name) in enumerate(opNames)]
# 		#remove the Identity operators since we don't need those
# 		operators = filter( (op) -> op.name != "I", operators)
# 		term = TERM(terms.prefactor, operators)
# 		push!(termsArray,term)
# 	end

# 	return Array{TERM}(termsArray)
# end




########################## new version #############################


"""

	SOP(name,site)

A struct for storing a single-site operator. This is an internal method which users shouldn't have to access.
"""
struct SOP
    name :: String
    site :: Int
end
Base.show(io::IO, sop::SOP) = print(io,  sop.name, "_", sop.site)


"""
	parse_term(op, term)

Internal function to parse the operators given in the `TERM` constructor into a Dictionary.
"""
function parse_term(
		operatorName :: String,
		firstSite :: Int
	)::Dict{Int, String}

	#get the names of all the operators
	#right now this assumes every character is separate operator
	#but eventually it should split differently based on the site type
	SOP_names = [string(ch) for ch in operatorName]
		
	operators = Dict( (firstSite+j-1) => sop_name for (j,sop_name) in enumerate(SOP_names))
	operators = filter( kv -> kv[2] != "I" && kv[2] != "1", operators)
	return operators

end

"""
	TERM(prefactor,operatorName,firstSite, [repeat=0])

	TERM(operatorName) = TERM(1, operatorName, 1, [repeat=1])
	TERM(prefactor,operatorName) = TERM(prefactor, operatorName,1, [repeat=1])
	TERM(operatorName,firstSite,[repeat=0]) = TERM(1, operatorName,1, [repeat=0])

A struct which stores one term in a Hamiltonian. A `TERM` is an operator whose action maps basis vectors to basis vectors (and not superpositions thereof). For spin-half systems, a term is any string of pauli operators, such as "Z_0 Z_1". The `TERM` struct can be constructed in several different ways fdepending on the situation.

# Arguments
- 'prefactor:: Number': the numerical prefactor for the term.
- 'operatorName:: String': a string of the allowed operators. For spin-half, the allowed operators are "1","I","X","Y","Z","+","-".
- 'firstSite:: Int': the site for the first Pauli operator in the string.
- 'repeat: Int':: the repeat period for the string of operators. If the `firstSite` is not specified, then the string repeats with period 1. If `firstSite` is specified, then the string does not repeat by default --- see examples below.

# Examples
```jldoctest
julia> TERM("ZZ")
TERM(1, "ZZ", 0, 1)

julia> TERM(3.0,"ZXXXZ",3)
TERM(3.0, "ZXXXZ", 3, 0)

julia> TERM("ZXZ",5,repeat=2)
TERM(1, "ZXZ", 5, 2)

julia> TERM(4.3,"ZXZ",5,repeat=2)
TERM(4.3, "ZXZ", 5, 2)

julia> TERM(-1.2, )
```
"""
struct TERM
	prefactor :: Number
	# operatorName :: String
	# firstSite :: Int
	operators :: Dict{Int,String}
	period :: Int

	function TERM(
		operatorName :: String
		)
		new(1,parse_term(operatorName,0),1)
	end

	function TERM(
		prefactor :: Number,
		operatorName :: String
		)
		new(prefactor,parse_term(operatorName,0),1)
	end

	function TERM(
		operatorName :: String,
		firstSite :: Int;
		repeat = 0 :: Int
		)
		@assert repeat >= 0 "The repeat period must be a positive integer."
		new(1,parse_term(operatorName,firstSite),repeat)
	end

	function TERM(
		prefactor :: Number,
		operatorName :: String,
		firstSite :: Int;
		repeat = 0 :: Int
		)
		@assert repeat >= 0 "The repeat period must be a positive integer."
		new(prefactor,parse_term(operatorName,firstSite),repeat)
	end

	function TERM(
		prefactor :: Number,
		operators :: Dict{Int, String};
		repeat = 0 :: Int)
		new(prefactor,operators,repeat)
	end
end



"""
	*(x,TERM)
	x*TERM

Defined scalar multiplication of `TERM`'s: the scalar acts on the prefactor.
"""
function Base.:*(x :: Number, t :: TERM)::TERM
    return TERM(x*t.prefactor,t.operators,repeat=t.period)
end

"""
	TERM_INTERNAL(prefactor, operator)

The internal representation of a term as a (generically complex) `prefactor` and an array of single-single operators called `operator`.
"""
struct TERM_INTERNAL
	prefactor :: Complex{Float64}
	operator :: Array{SOP}
end
function Base.show(io::IO, term::TERM_INTERNAL)
	print(io,  term.prefactor, "*")
	for sop in term.operator
		print(io,sop," ")
	end
end

"""
	*(x,TERM_INTERNAL)
	x*TERM

Defined scalar multiplication of `TERM`'s: the scalar acts on the prefactor.
"""
function Base.:*(x :: Number, t :: TERM_INTERNAL)::TERM_INTERNAL
    return TERM_INTERNAL(x*t.prefactor,t.operators)
end


# ABSTRACT_OP(L, sites, name, pbc, terms)
"""

	ABSTRACT_OP(L; sites = "spin half", name = "abstract operator", pbc = true)
	ABSTRACT_OP(L, operatorName, site; sites = "spin half", name = "abstract operator", pbc = true)
	ABSTRACT_OP(L, TERM; sites = "spin half", name = "abstract operator", pbc = true)

The `ABSTRACT_OP` struct represents an operator as a string of operator names on sites, along with their numerical prefactors. This also encodes some details of the Hilbert space, including the number of sites `L`, the type of site in use `sites` (e.g. spin half), the name for the operator `name`, and `pbc` which is true when periodic boundary conditions are employed.

`ABSTRACT_OP`'s constitute a vector space and can be added and scalar-multiplied by (generically) complex numbers. 

# Examples

There are several different types of constructors available for `ABSTRACT_OP`, enabling the quick definition of "shells" of operators which can be used to define complex Hamiltonians, or quick constructors for simple observables.

```jldoctest
julia> ABSTRACT_OP(10)
ABSTRACT_OP[name: abstract operator, L: 10, type: spin half, pbc: true, #terms: 0]

julia> ABSTRACT_OP(10,"X",4)
ABSTRACT_OP[name: abstract operator, L: 10, type: spin half, pbc: true, #terms: 1]
1.0 + 0.0im*X_4

julia> ABSTRACT_OP(10,TERM("Y",5))
ABSTRACT_OP[name: abstract operator, L: 10, type: spin half, pbc: true, #terms: 1]
1.0 + 0.0im*Y_5

julia> ABSTRACT_OP(10,TERM("Z",7); pbc=false)
ABSTRACT_OP[name: abstract operator, L: 10, type: spin half, pbc: false, #terms: 1]
1.0 + 0.0im*Z_7

julia> ABSTRACT_OP(10,4.3TERM("Z",7); name="S_z^7", pbc=false)
ABSTRACT_OP[name: S_z^7, L: 10, type: spin half, pbc: false, #terms: 1]
4.3 + 0.0im*Z_7

julia> H = ABSTRACT_OP(4; name="Ising Model", pbc=true)
julia> H += TERM("ZZ")
julia> H += TERM("X")
ABSTRACT_OP[name: Ising Model, L: 4, type: spin half, pbc: true, #terms: 8]
1.0 + 0.0im*Z_0 Z_1
1.0 + 0.0im*Z_1 Z_2
1.0 + 0.0im*Z_2 Z_3
1.0 + 0.0im*Z_3 Z_0
1.0 + 0.0im*X_0
1.0 + 0.0im*X_1
1.0 + 0.0im*X_2
1.0 + 0.0im*X_3

julia> order_parameter = ABSTRACT_OP(4,0.25*TERM("ZZ"))
ABSTRACT_OP[name: abstract operator, L: 4, type: spin half, pbc: true, #terms: 4]
0.25 + 0.0im*Z_0 Z_1
0.25 + 0.0im*Z_1 Z_2
0.25 + 0.0im*Z_2 Z_3
0.25 + 0.0im*Z_3 Z_0
```

See also: [`TERM`](@ref).
"""
struct ABSTRACT_OP
	L :: Int
	sites :: String
    name :: String
    pbc :: Bool
    terms :: Array{TERM_INTERNAL}

    function ABSTRACT_OP(
    	L:: Int;
    	sites = "spin half" :: String,
    	name = "abstract operator" :: String,
    	pbc = true :: Bool
    	)
    	new(L,sites,name,pbc,[])
	end

	function ABSTRACT_OP(
    	L:: Int,
    	operatorName :: String,
    	site :: Int;
    	sites_type = "spin half" :: String,
    	name = "abstract operator" :: String,
    	pbc = true :: Bool
    	)
    	new(L,sites_type,name,pbc,[]) + TERM(operatorName,site)
	end

	function ABSTRACT_OP(
    	L:: Int,
    	term :: TERM; 
    	sites_type = "spin half" :: String,
    	name = "abstract operator" :: String,
    	pbc = true :: Bool
    	)
    	new(L,sites_type,name,pbc,[]) + term
	end

	function ABSTRACT_OP(
    	L:: Int,
    	sites :: String,
    	name :: String,
    	pbc :: Bool,
    	terms :: Array{TERM_INTERNAL}
    	)
    	new(L,sites,name,pbc,terms)
	end
end

function Base.show(io::IO, op::ABSTRACT_OP)
	print(io, "ABSTRACT_OP[name: \"$(op.name)\", L: $(op.L), type: $(op.sites), pbc: $(op.pbc), #terms: $(length(op.terms))]")
	for t in op.terms
		print(io,"\n", t)
	end
end


#ABSTRACT_OP's should form a vector space with
#addition and scalar multiplication

"""

	O1 + O2

Adds two abstract operators. Operators must have the same length, site type, and boundary conditions.
"""
function Base.:+(O1 :: ABSTRACT_OP, O2 :: ABSTRACT_OP)::ABSTRACT_OP
	@assert O1.L == O2.L "Operators must have the same length to be added."
	@assert O1.sites == O2.sites "Operators must have the same site type to be added."
	@assert O1.pbc == O2.pbc "Operators must have the same boundary conditions to be added."

	newName = string(O1.name," + ",O2.name)
	newTerms = [O1.terms;O2.terms]
	return ABSTRACT_OP(O1.L,O1.sites,newName,O1.pbc,newTerms)
end

"""
	x*Op

Scalar multiplication of operators.
"""
function Base.:*(lambda :: Complex{Float64}, O :: ABSTRACT_OP)::ABSTRACT_OP
	newTerms = [TERM_INTERNAL(lambda*t.prefactor,t.operator) for t in O.terms]
	return ABSTRACT_OP(O.L,O.sites,O.name,O.pbc,newTerms)
end


"""
	op + term

Adds a new `TERM` to an `ABSTRACT_OP`. The `TERM` must fit inside the number of sites for the operator.

# Examples
```jldoctest

julia> H = ABSTRACT_OP(4; name="Ising Model", pbc=true)
julia> H += TERM("ZZ")
julia> H += TERM("X")
ABSTRACT_OP[name: Ising Model, L: 4, type: spin half, pbc: true, #terms: 8]
1.0 + 0.0im*Z_0 Z_1
1.0 + 0.0im*Z_1 Z_2
1.0 + 0.0im*Z_2 Z_3
1.0 + 0.0im*Z_3 Z_0
1.0 + 0.0im*X_0
1.0 + 0.0im*X_1
1.0 + 0.0im*X_2
1.0 + 0.0im*X_3
```


"""
function Base.:+(op :: ABSTRACT_OP, term :: TERM)::ABSTRACT_OP
	#terminate early if the term is zero to floating-point precision
	if abs(term.prefactor) <= 10e-15
		return op
	else	 
		newTerms = [op.terms; parse_TERM_to_internal(op,term)]
		return ABSTRACT_OP(op.L,op.sites,op.name,op.pbc,newTerms)
	end
end

"""
	parse_term(op, term)

Internal function to parse a `TERM` into an array of `TERM_INTERNAL`s so they can be added onto an operator.
"""
function parse_TERM_to_internal(
		op :: ABSTRACT_OP,
		term :: TERM
	)::Array{TERM_INTERNAL}

	#harcoded spin-half for now, but eventually this could be expanded to other types of sites
	allowed_operatorNames_spin_half = ["1","I","X","Y","Z","+","-"]
	# println(allowed_operatorNames_spin_half)

	#check all the names are valid
	#change to a hashmap?
	for (site,sop_name) in term.operators
		@assert in(sop_name, allowed_operatorNames_spin_half) "Unknown operator name \"$(sop_name)\" in , $(term)"
	end 
	# println("single operator names:", SOP_names)

	first_site = minimum(term.operators)[1]
	last_site = maximum(term.operators)[1]

	#check that at least one operator fits
	@assert first_site >= 0 && last_site < op.L "Operator $(term) does not fit. Site range: 0 to $(op.L-1)"


	#make the array of starting sites
	termEnd = op.pbc ?  op.L-1-first_site : (op.L)-last_site-1
	# println("termEnd:", termEnd)
	offsets = term.period == 0 ? [0] : collect(0 : term.period : termEnd)
	# println("offsets:",offsets)

	#actually make the terms now
	termsArray = Array{TERM_INTERNAL}(undef,length(offsets))

	#make the array of SOP's for each starting site
	for (n,offset) in enumerate(offsets)
		operators = [SOP(sop_name,(site+offset)%op.L) for (site,sop_name) in term.operators]
		operators = sort(operators, by = sop -> sop.site)

		#repackage as an internal term struct
		termsArray[n] = TERM_INTERNAL(term.prefactor, operators)
	end

	return termsArray
end