#####################################################################
#	This file implements the idea of "abstract operators":
#	abstract representations of operators as strings of Pauli matrices
#	or other symbolic representations, which can easily be turned into 
#	real matrices for a given Hilbert space.
#
#	The main data structures are ABSTRACT_OP, which represents such an operator,
#	and TERM, which gives a human-readable interface to build them
#	in the same way that spin-chains are written on paper.
#####################################################################


"""

	SOP(name,site)

A struct for storing a single-site operator. This is an internal method which users shouldn't have to access.
"""
struct SOP
    name :: String
    site :: Int
end
Base.show(io::IO, sop::SOP) = print(io,  sop.name, "_", sop.site)

# """
# [Update doc]

# A function for ordering single-site operators.
# """
# function isless(SOPA :: SOP, SOPB :: SOP)
# 	return SOPA.site < SOPB.site
# end


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
	TERM(prefactor,operatorName, firstSite, [repeat=0])

	TERM(operatorName) = TERM(1, operatorName, 1, [repeat=1])
	TERM(prefactor,operatorName) = TERM(prefactor, operatorName,1, [repeat=1])
	TERM(operatorName,firstSite,[repeat=0]) = TERM(1, operatorName,1, [repeat=0])
	TERM(prefactor, operatorDictionary; repeat = 0)

A struct which stores one term in a Hamiltonian. A `TERM` is an operator whose action maps basis vectors to basis vectors (and not superpositions thereof). For spin-half systems, a term is any string of pauli operators, such as ``Z_0 Z_1``. The `TERM` struct can be constructed in several different ways fdepending on the situation.

# Arguments
- 'prefactor:: Number': the numerical prefactor for the term.
- 'operatorName:: String': a string of the allowed operators. For spin-half, the allowed operators are "1", "I", "X", "Y", "Z", "+", and "-".
- 'firstSite:: Int': the site for the first Pauli operator in the string.
- 'repeat: Int':: the repeat period for the string of operators. If the `firstSite` is not specified, then the string repeats with period 1. If `firstSite` is specified, then the string does not repeat by default --- see examples below.

# Examples
```jldoctest
julia> TERM("ZZ")
TERM(1, Dict(0=>"Z",1=>"Z"), 1)

julia> TERM(3.0,"ZXXXZ",3)
TERM(3.0, Dict(7=>"Z",4=>"X",3=>"Z",5=>"X",6=>"X"), 0)

julia> TERM("ZXZ",5; repeat=2)
TERM(1, Dict(7=>"Z",5=>"Z",6=>"X"), 2)

julia> TERM(4.3,"ZXZ",5,repeat=2)
TERM(4.3, Dict(7=>"Z",5=>"Z",6=>"X"), 2)

julia> TERM(0.5,Dict(1=>"X", 7=>"X"))
TERM(0.5, Dict(7=>"X",1=>"X"), 0)
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
	x * TERM
	*(x,TERM)

Scalar multiplication of `TERM`'s: the scalar 'x' acts on the prefactor to 'TERM'.

# Examples

```jldoctest
julia> t1 = TERM(4,"ZXZ",5,repeat=2)
TERM(4, Dict(7=>"Z",5=>"Z",6=>"X"), 2)
julia> 3im*t1
TERM(0 + 12im, Dict(7=>"Z",5=>"Z",6=>"X"), 2)

```
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

#####
##### for now, removed references to the `sites` option since it's all spin-1/2
#####
# 	ABSTRACT_OP(L; sites = "spin half", name = "abstract operator", pbc = true)
# 	ABSTRACT_OP(L, operatorName, site; sites = "spin half", name = "abstract operator", pbc = true)
# 	ABSTRACT_OP(L, TERM; sites = "spin half", name = "abstract operator", pbc = true)

# The `ABSTRACT_OP` struct represents an operator as a string of operator names on sites, along with their numerical prefactors. This also encodes some details of the Hilbert space, including the number of sites `L`, the type of site in use `sites` (e.g. spin half), the name for the operator `name`, and `pbc` which is true when periodic boundary conditions are employed.
# ABSTRACT_OP(L, sites, name, pbc, terms)
"""

	ABSTRACT_OP(L; name = "abstract operator", pbc = true)
	ABSTRACT_OP(L, operatorName, site; name = "abstract operator", pbc = true)
	ABSTRACT_OP(L, TERM; name = "abstract operator", pbc = true)

The `ABSTRACT_OP` struct represents an operator as a string of operator names on sites, along with their numerical prefactors. This also encodes some details of the Hilbert space, including the number of sites `L`, the name for the operator `name`, and `pbc` which is true when periodic boundary conditions are employed.

`ABSTRACT_OP`'s constitute a vector space and can be added and scalar-multiplied by (generically) complex numbers. 

# Examples

There are several different types of constructors available for `ABSTRACT_OP`, enabling the quick definition of "shells" of operators which can be used to define complex Hamiltonians, or quick constructors for simple observables.

```jldoctest
julia> ABSTRACT_OP(10)
ABSTRACT_OP[name: "abstract operator", L: 10, type: spin half, pbc: true, #terms: 0]

julia> ABSTRACT_OP(10,"X",4)
ABSTRACT_OP[name: "abstract operator", L: 10, type: spin half, pbc: true, #terms: 1]
1.0 + 0.0im*X_4

julia> ABSTRACT_OP(10,TERM("Y",5))
ABSTRACT_OP[name: "abstract operator", L: 10, type: spin half, pbc: true, #terms: 1]
1.0 + 0.0im*Y_5

julia> ABSTRACT_OP(10,TERM("Z",7); pbc=false)
ABSTRACT_OP[name: "abstract operator", L: 10, type: spin half, pbc: false, #terms: 1]
1.0 + 0.0im*Z_7

julia> ABSTRACT_OP(10,4.3TERM("Z",7); name="S_z^7", pbc=false)
ABSTRACT_OP[name: "S_z^7", L: 10, type: spin half, pbc: false, #terms: 1]
4.3 + 0.0im*Z_7

julia> ABSTRACT_OP(4; name="Ising Model", pbc=true) + TERM("ZZ") + TERM("X")
ABSTRACT_OP[name: "Ising Model", L: 4, type: spin half, pbc: true, #terms: 8]
1.0 + 0.0im*Z_0 Z_1
1.0 + 0.0im*Z_1 Z_2
1.0 + 0.0im*Z_2 Z_3
1.0 + 0.0im*Z_0 Z_3
1.0 + 0.0im*X_0
1.0 + 0.0im*X_1
1.0 + 0.0im*X_2
1.0 + 0.0im*X_3

julia> order_parameter = ABSTRACT_OP(4,0.25*TERM("ZZ"))
ABSTRACT_OP[name: "abstract operator", L: 4, type: spin half, pbc: true, #terms: 4]
0.25 + 0.0im*Z_0 Z_1
0.25 + 0.0im*Z_1 Z_2
0.25 + 0.0im*Z_2 Z_3
0.25 + 0.0im*Z_0 Z_3
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
	+(O1,O2)

Adds two `ABSTRACT_OP`s. Operators must have the same length, site type, and boundary conditions.

# Examples
```jldoctest
julia> 	OP1 = ABSTRACT_OP(10,"X",4); OP2 = ABSTRACT_OP(10,"Z",5);

julia> OP1+OP2
ABSTRACT_OP[name: "abstract operator + abstract operator", L: 10, type: spin half, pbc: true, #terms: 2]
1.0 + 0.0im*X_4
1.0 + 0.0im*Z_5
```
"""
function Base.:+(O1 :: ABSTRACT_OP, O2 :: ABSTRACT_OP)::ABSTRACT_OP
	@assert O1.L == O2.L "Operators must have the same length to be added."
	@assert O1.sites == O2.sites "Operators must have the same site type to be added."
	@assert O1.pbc == O2.pbc "Operators must have the same boundary conditions to be added."

	newName = string(O1.name," + ",O2.name)
	newTerms = [O1.terms;O2.terms]
	return ABSTRACT_OP(O1.L,O1.sites,newName,O1.pbc,newTerms)
end

function Base.:-(O1 :: ABSTRACT_OP, O2 :: ABSTRACT_OP)::ABSTRACT_OP
	O1 + (-O2)
end

"""
	x * Op
	*(x, Op)

Scalar multiplication of `ABSTRACT_OP`s.
# Examples
```jldoctest
julia> 	Op = ABSTRACT_OP(10,TERM(2,"X",4))
ABSTRACT_OP[name: "abstract operator", L: 10, type: spin half, pbc: true, #terms: 1]
2.0 + 0.0im*X_4

julia> 3*Op
ABSTRACT_OP[name: "abstract operator", L: 10, type: spin half, pbc: true, #terms: 1]
6.0 + 0.0im*X_4

```
"""
function Base.:*(lambda :: Number, O :: ABSTRACT_OP)::ABSTRACT_OP
	newTerms = [TERM_INTERNAL(lambda*t.prefactor,t.operator) for t in O.terms]
	return ABSTRACT_OP(O.L,O.sites,O.name,O.pbc,newTerms)
end


"""
	Op + term
	+(Op, term)

Adds a new `TERM` to an `ABSTRACT_OP`. The `TERM` must fit inside the number of sites for the operator.

# Examples
```jldoctest
julia> ABSTRACT_OP(4; name="Ising Model", pbc=true) + TERM("ZZ") + TERM("X")
ABSTRACT_OP[name: "Ising Model", L: 4, type: spin half, pbc: true, #terms: 8]
1.0 + 0.0im*Z_0 Z_1
1.0 + 0.0im*Z_1 Z_2
1.0 + 0.0im*Z_2 Z_3
1.0 + 0.0im*Z_0 Z_3
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

Base.:-(op :: ABSTRACT_OP, term :: TERM)::ABSTRACT_OP = op + (-term)

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


"""
[Update doc]

A function for consolidating multiple coefficients for the same operator.
"""
function consolidate_coefficients(Op :: ABSTRACT_OP)

	termsDict = Dict{Array{SOP},Complex{Float64}}()

	for t in Op.terms
		sorted = sort(t.operator,by= x-> x.site)
		if !haskey(termsDict,sorted)
			termsDict[sorted] = t.prefactor
		else
			termsDict[sorted] += t.prefactor
		end
	end

	termsArray = [TERM_INTERNAL(factor,term) for (term,factor) in termsDict]

	return ABSTRACT_OP(Op.L,Op.sites,Op.name,Op.pbc,termsArray)
end