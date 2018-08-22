#################### Data Structures ############################

"""
	BASISVECTOR(conj_class, phase_factor)

A struct that links a basis vector to its conjugacy class `conj_class' and gives the phase factor 'phase_factor' between them due to the representation.

"""
struct BASISVECTOR
    conj_class :: UInt64
    phase_factor :: ComplexF64
end
Base.show(io::IO, bv::BASISVECTOR) = print(io, "BV[",bv.conj_class, ", ", bv.phase_factor, "]")

"""
	CONJCLASS(index,norm)

A struct that describes the conjugacy class of a basis element. With symmetries, each conjugacy class is associated to a basis vector. The `index` is the index of the associated basis vector and the `norm` is the *norm squared* of that conjugacy class, so it can be properly normalized.
"""
struct CONJCLASS
    index :: Int64
    norm :: Int
end
Base.show(io::IO, cc::CONJCLASS) = print(io, "cc[ind:", cc.index,", norm:", cc.norm, "]")

"""
	ORBIT(norm, representative, elements)

A struct which stores the orbit of a single element under an Abelian group action.
"""
struct ORBIT
	norm :: Int
	representative :: UInt64
	elements :: Dict{UInt64,ComplexF64}
end
Base.show(io::IO,orb::ORBIT) = print(io,
"ORBIT[",
"\n\tnorm: ", orb.norm,
"\n\trep: ", orb.representative,
"\n\telements: \n",
[ "$(string(k,base=2)) => $(zchop(v))" for (k,v) in orb.elements],
"]"
	)

##################### Symmetry operations  #########################

"""
	measure_N_up(s, L)

Measures the number of spin up's for a state 's' for 'L' spins.
Returns the Sz sector of s, i.e. the number of spin's up.
"""
function measure_N_up(s::UInt64, L::Int)
    #s - state
    #L - Length of the spin chain, must be even
    #returns - the Sz sector of s, i.e. the number of spins up
    a = Int(sum([(s & (1 << pl)) >> pl for pl = 0:L-1]))
    return a
end

"""
	measure_N_up_A(s, L)

Measures the number of spin up's for a state 's' for 'L' spins in the A sector,
for an ABABAB sublattice structure.
Returns the SzA sector of s.
"""
function measure_N_up_A(s :: UInt64, L :: Int)
    return Int(sum([(s & (1 << pl)) >> pl for pl = 0:2:L-1]))
end

"""
	measure_N_up_B(s, L)

Measures the number of spin up's for a state 's' for 'L' spins in the B sector,
for an ABABAB sublattice structure.
Returns the SzB sector of s.
"""
function measure_N_up_B(s :: UInt64, L :: Int)
    return Int(sum([(s & (1 << pl)) >> pl for pl = 1:2:L-1]))
end

################ Functions to apply symmetry ###############

"""

	make_translation_function(L,a,K)
	
Returns a function which, given a state x and phase factor pf,
gives the translated state T.x and new phase factor.

#Arguments
* 'L :: Int': the number of sites
* 'a :: Int': the number of sites per unit cell, must divide L
* 'K :: Int': the translation symmetry sector, 0 <= K <= L
"""
function make_translation_function(L :: Int, a :: Int,  K :: Int)
	G_k_size = div(L,a)
	omega :: ComplexF64 = zchop(exp(- (2 * pi * im * K)/G_k_size))
	c1 = UInt64(2^L-1)

	return 	function (x :: UInt64, pf :: ComplexF64)
				shift :: UInt64 = x << a
				gx :: UInt64 = (shift & c1) | ((shift & ~c1) >> L)
				pf *= omega

				return (gx, pf)
			end
end

"""

	make_Z2B_function(L,Z2B)

Returns a function which computes a version of a state flipped on the B sublattice. This only makes sense when unitCellSize = 2, of course.
#Arguments
* 'L :: Int': the number of sites
* 'Z2B :: Int': the Z2B sector, in {-1,1}
"""
function make_Z2B_function(
	L :: Int,
	Z2B :: Int
	)
	flipper :: UInt64 = UInt64(sum([1 << (k+1) for k in 0:2:L-2]))

	if Z2B == 1
		return 	function (x :: UInt64, pf :: ComplexF64)
					gx :: UInt64 = xor(x, flipper)
					return (gx,pf)
				end
	else 
		return 	function (x :: UInt64, pf :: ComplexF64)
							gx :: UInt64 = xor(x,flipper)
							pf *= -1
							return (gx,pf)
				end
	end
end


"""

	make_Z2A_function(L,Z2A)

Returns a function which computes a version of a state flipped on the B sublattice. This only makes sense when unitCellSize = 2, of course.
#Arguments
* 'L :: Int': the number of sites
* 'Z2A :: Int': the Z2A sector, in {-1,1}
"""
function make_Z2A_function(
	L :: Int,
	Z2A :: Int
	)

	flipper :: UInt64 = UInt64(sum([1 << k for k in 0:2:L-1]))

	if Z2A == 1
		return 	function (x :: UInt64, pf :: ComplexF64)
					gx :: UInt64 = xor(x, flipper)
					return (gx,pf)
				end
	else 
		return 	function (x :: UInt64, pf :: ComplexF64)
							gx :: UInt64 = xor(x,flipper)
							pf *= -1
							return (gx,pf)
						end
	end
end


"""

	make_spin_flip_function(L,Z2)

Returns a function which computes the spin flip of a state
#Arguments
* 'L :: Int': the number of sites
* 'Z2 :: Int': the Z2 sector, in {-1,1}
"""
function make_spin_flip_function(
	L :: Int,
	Z2 :: Int
	)
	
	flipper :: UInt64 = UInt64(sum([1 << (k) for k in 0:L-1]))
	# println("flipper: ", string(flipper,base=2))

	if Z2 == 1
		return 	function (x :: UInt64, pf :: ComplexF64)
					gx :: UInt64 = xor(x, flipper)

					return (gx,pf)
				end
	else 
		return 	function (x :: UInt64, pf :: ComplexF64)
							gx :: UInt64 = xor(x,flipper)
							pf *= -1
							return (gx,pf)
						end
	end
end

"""

	make_Inversion_function(L,Inv,a)

Returns a function which computes the inversion (flip in the x-direction) of a state
#Arguments
* 'L :: Int': the number of sites
* 'Inv :: Int': the Inversion sector, in {-1,1}
* 'a :: Int': the size of the unit cell, must divide L
"""
function make_Inversion_function(
	L :: Int,
	Inv :: Int,
	a :: Int
	)
	
	cell = UInt64((2^a) -1)

	#this seems fairly inefficient, but I can't think of a better way right now
	if Inv == 1
		return 	function (x :: UInt64, pf :: ComplexF64)
					gx = UInt64(0)
					shift_back = (L-a)
					for k in 0:a:L-1
						k_cell = cell << k #kth cell mask
						x2 = x & k_cell	 	#extract kth cell only
						x3 = x2 << shift_back #push it back
						shift_back -= 2a	#update shift
						gx = gx | x3		#update vector
						#println(k, ", x2: ", x2, " x3: ", x3, " gx: ", gx)
					end
					return (gx,pf)
				end
	else 
		return 	function (x :: UInt64, pf :: ComplexF64)
					gx = UInt64(0)
					shift_back = (L-a)
					for k in 0:a:L-1
						k_cell = cell << k #kth cell mask
						x2 = x & k_cell	 	#extract kth cell only
						x3 = x2 << shift_back #push it back
						shift_back -= 2a	#update shift
						gx = gx | x3		#update vector
						#println(k, ", x2: ", x2, " x3: ", x3, " gx: ", gx)
					end
					pf *= -1
					return (gx,pf)
				end
	end
end



############ Make a "universal" basis ##################


####symmetries are broken into one of two types:
# 1. symmetries that restrict which states are valid,
# 	 such as magnetization sectors
# 2. and "proper" symmetries with associated phase factors
# The code filters the first out separately, since they're much faster



"""
	get_valid_state_function(L, unitCellSize,, symmetries)


Returns a function  `states -> Bool` that determines if a state is valid, i.e. satisfies the constraint imposed by the symmetry sector.

# Arguments
- 'L :: Int': length of the chain
- 'unitCellSize:: Int' - size of the unit cell, must divide L
- 'symmetries :: Dict{String,Int}': a dictionary of the symmetries and their sectors, e.g ["Sz"=> 2,"K"=>2]

"""
function get_valid_state_function(L :: Int, unitCellSize :: Int, symmetries :: Dict{String,Int})
    #L - length
    #unitCellSize - size of the unit cell
    #symmetries - Array of the symmetries, e.g [("Sz", 2),("K",-2)]
    #returns - a function that determines if a state is valid

    if haskey(symmetries,"Sz")
		N_up = symmetries["Sz"]
    	return x -> measure_N_up(x,L) == N_up
    elseif haskey(symmetries,"SzA") && unitCellSize == 2
    	N_up_A = symmetries["SzA"]
    	return x -> measure_N_up_A(x,L) == N_up_A
    elseif haskey(symmetries,"SzB") && unitCellSize == 2
    	N_up_B = symmetries["SzB"]
    	return x -> measure_N_up_B(x,L) == N_up_B
   	else
   		return x -> true
    end 
end


"""

	make_easy_orbit_function(L,unitCellSize,symmetries)	

Given the symmetry information, returns a function which computes the orbit of an individual element. Works for valid basis elements only.

# Arguments
- 'L :: Int': length of the chain
- 'unitCellSize:: Int' - size of the unit cell, must divide L
- 'symmetries :: Dict{String,Int}': a dictionary of the symmetries and their sectors, e.g ["Sz"=> 2,"K"=>2]

"""
function make_easy_orbit_function(
	L :: Int,
	unitCellSize :: Int,
	symmetries :: Dict{String,Int}
	)

	symmetries_fcns = Vector()
	G_size = 1

	for (name,sector) in symmetries
		if name == "PA"
			sym_fcn = make_Z2A_function(L,sector)
			sub_gp_size = 2
		elseif name == "PB"
			sym_fcn = make_Z2B_function(L,sector)
			sub_gp_size = 2
		elseif name == "I"
			sym_fcn = make_Inversion_function(L,sector,unitCellSize)
			sub_gp_size = 2
		elseif name == "P"
			sym_fcn = make_spin_flip_function(L,sector)
			sub_gp_size = 2
		elseif name == "Tr"
			@assert L % unitCellSize == 0 "Unit cell size $(unitCellSize) does not divide the number of spins $(L)."
			sym_fcn = make_translation_function(L,unitCellSize,sector)
			sub_gp_size = div(L,unitCellSize)
		end
		G_size *= sub_gp_size
		push!(symmetries_fcns, (sub_gp_size, sym_fcn))
	end			

	# println("symmetries are:")
	# println(symmetries_fcns)

	return function (x :: UInt64)
		cc_x = x
		norm_x :: ComplexF64 = 1.0
		phase_factor = Complex(1.0)
		O_x =  Dict{UInt64, ComplexF64}(x => phase_factor)

		for (subgp_size, sym_fcn) in symmetries_fcns
			norm_factor = norm_x # how many multiples do we have already
			O_x_new = Dict(O_x)
			if subgp_size == 2
				#println("Debug: applying non-translation symmetry.")
				for (y,pf) in O_x
					(gy, pf_new) = sym_fcn(y,pf)
					if gy == x
						norm_x += pf_new*norm_factor
					elseif !haskey(O_x_new,gy)
						O_x_new[gy] = pf_new
						if gy < cc_x
							cc_x = gy
						end
					end
				end
			else#this must be translation, at this point
				#println("Debug: applying translation.")
				for (y,pf) in O_x
					for i in 2:subgp_size
						(gy, pf_new) = sym_fcn(y,pf)
						if gy == x
							norm_x += pf_new*norm_factor
						elseif !haskey(O_x_new,gy)
							O_x_new[gy] = pf_new
							if gy < cc_x
								cc_x = gy
							end
						end
						y = gy
						pf = pf_new
					end
				end
			end

			if abs(norm_x) <= 10e-14#terminate early if the norm is already zero
				return Nullable{ORBIT}()
			end
			O_x = O_x_new
		end
		orb = ORBIT(div(G_size,length(O_x)), cc_x, O_x)
		return Nullable(orb)
	end
end
# ###eventually it would be good to come back to this and optimize for efficiency,
# #e.g. make the inner loops into their own funcitons
# """
# 	make_basis(L; [unitCellSize=a, syms=Symmetries_Dict, constraint=constraint_function])

# Produces a basis from its symmetry information and a function that constrains the Hilbert space.
# # Arguments
# - 'L :: Int': the number of sites
# - 'unitCellSize :: Int': the number of sites per unit cell
# - 'syms :: Dict{String,Int}': a dictionary of symmetries with their sector names, e.g. the translation sector 3 is "K" => 3
# - 'constraint :: Function': a function (state, L) -> bool to determine if a given state satisfies a constraint
# """

####Type for storing bases, together with its constructor

"""
		BASIS(L;
		a :: Int64 = 1,
		syms :: Dict{String,Int} = (),
		constraint :: Function = identity
		)

A struct that stores a basis for a Hilbert space. A `BASIS` is constructed by specifiying its attributes.
	- The number of spins `L`
	- The number of spins per unit cell `a`, whose default is 1.
	- Any symmetries of the basis, specified as a dictionary `syms`. The default is no symmetries.
	- A constraint function `state :: UInt64 -> Bool` where only states where the constraint evaluates to `true` are admitted to the Hilbert space.

Symmetries are specified by giving their name and sector. For instance, `"tr" => 3` represents the translation symmetry in sector 3, where one application of the symmetry moves the state around by 3 basis vectors. The complete list of implemented symmetries is below.

| Symmetry   | Full Name                            | Values                    |
|:----------:|--------------------------------------|:-------------------------:|
| Tr         | Translation                          | ``1 \\le Tr \\le L/a``    |
| P          | Parity (Spin flip)                   | ``\\{-1,1\\}``            |
| I          | Inversion                            | ``\\{-1,1\\}``            |
| Sz         | Total ``S^z``                        | ``-L \\le Sz \\le L``     |
| PA         | P for even spins                     | ``\\{-1,1\\}``            |
| PB         | P for odd spins                      | ``\\{-1,1\\}``            |
| SzA        | ``\\sum_{i\\;\\mathrm{even}} S^z_i`` | ``-L/2 \\le SzA \\le L/2``|
| SzB        | ``\\sum_{i\\;\\mathrm{odd}} S^z_i``  | ``-L/2 \\le SzB \\le L/2``|

Using some symmetries requires extra conditions. Explicitly, translation symmetry requires periodic boundary conditions, and even/odd parity and spin sectors require a even number of total sites. However, multiple symmetries can be used at the same time.

!!! warning "Index Convention"

    Spins are zero-indexed, so the left-most site is spin zero and is in the "A" sublattice.

Users can also supply arbitrary constraint functions `UInt64 -> Bool`, which specify states that are allowed in the Hilbert space. 

!!! note

    A `BASIS` can require a very large amount of memory to store. Internally, it is composed of one hashmap from states in the full Hilbert space to a representative of their conjugacy classes, and another from representatives of conjugacy classes to indices in the small Hilbert space.
    Altogether, this requires 2^L space. For truly large sizes, one can --- in principle --- generate these hashmaps as needed, but this is not implemented yet here.

# Examples
```jldoctest
julia> BASIS(8)
BASIS[L: 8, a: 1, dim: 256, symmetries: none]

julia> BASIS(8; a=2)
BASIS[L: 8, a: 2, dim: 256, symmetries: none]

julia> BASIS(8; syms = Dict("Tr"=>1))
BASIS[L: 8, a: 1, dim: 30, symmetries: "Tr" => 1]

julia> BASIS(8; a=2, syms = Dict("Tr"=>1))
BASIS[L: 8, a: 2, dim: 60, symmetries: "Tr" => 1]

julia> BASIS(8; syms = Dict("Sz" => 3, "Tr" => 2))
BASIS[L: 8, a: 1, dim: 7, symmetries: "Sz" => 3, "Tr" => 2]

julia> cons_fcn = x -> x < 17;

julia> BASIS(8; constraint = cons_fcn)
BASIS[L: 8, a: 1, dim: 17, symmetries: "constrained" => 1]
```
"""
struct BASIS
    #type for storing a basis with various quantum numbers
    L :: Int #number of sites
    a :: Int
    dim :: Int
    get_conj_class :: Dict{UInt64,BASISVECTOR} #gets basis vector for any state
    conj_classes :: Dict{UInt64,CONJCLASS}    #gets index and norm of a conjugacy class
    q_numbers :: Dict{String,Int} #list of quantum numbers for the state
    #e.g. q_numbers = [("K",1),("Sz",3)]

    #expose default constuctor
    function BASIS(
    	L :: Int,
    	a :: Int,
    	dim :: Int,
	    get_conj_class :: Dict{UInt64,BASISVECTOR}, #gets basis vector for any state
	    conj_classes :: Dict{UInt64,CONJCLASS},    #gets index and norm of a conjugacy class
	    q_numbers :: Dict{String,Int}
	    )
    	new(
		    L,
	    	a,
	    	dim,
	    	get_conj_class,
	    	conj_classes,
	    	q_numbers)
    end

	function BASIS(
			L :: Int;
			a :: Int64 = 1,
			syms :: Dict{String,Int} =  Dict{String,Int}(),
			constraint :: Function = identity
		)
		make_basis(L,unitCellSize=a,syms=syms,constraint=constraint)
	end
end

function Base.show(
	io::IO,
	b::BASIS)
	print(io, "BASIS[L: $(b.L), a: $(b.a), dim: $(b.dim), symmetries: ")
	if b.q_numbers != Dict{String,Int64}()
		for (n,(k,v)) in enumerate(b.q_numbers)
			if n != 1
				print(io,", ")
			end
			print(io,"\"$(k)\" => $(v)")
		end
	else
		print(io,"none")
	end
	print(io,"]")
end



###eventually it would be good to come back to this and optimize for efficiency,
#e.g. make the inner loops into their own funcitons
"""
	make_basis(L; [unitCellSize=a, syms=Symmetries_Dict, constraint=constraint_function])

Internal function to produce a basis from its size, constraints, and symmetry information.

# Arguments
- 'L :: Int': the number of sites
- 'unitCellSize :: Int': the number of sites per unit cell
- 'syms :: Dict{String,Int}': a dictionary of symmetries with their sector names, e.g. the translation sector 3 is "K" => 3
- 'constraint :: Function': a function (state, L) -> bool to determine if a given state satisfies a constraint
"""
function make_basis(
	L :: Int;  						#number of lattice sites
	unitCellSize :: Int = 1, 			#self-expanatory
	# validity_syms :: Dict{String,Int} =  Dict{String,Int}(),#symmetries that restrict which states are valid
	# #name => sector, e.g. translation sector 3 is "K" => 3
	syms :: Dict{String,Int} =  Dict{String,Int}(),
	constraint :: Function = identity
	):: BASIS

	@assert L % unitCellSize == 0 "The unit cell size, $(unitCellSize), does not divide the number of spins, $(L)"

	if unitCellSize == 1 && syms == Dict{String,Int}() && constraint === identity
		return BASIS(L,
			1, 
			2^L,
			Dict{UInt64,BASISVECTOR}(),
			Dict{UInt64,CONJCLASS}(),
			Dict{String,Int}())
	end

	#separate out the symemtries that simply restrict the Hilbert space
	#hardcoded names for now, as only a few are implemented
	validity_syms_names = ["Sz","SzA","SzB"]
	symmetries_names = ["Tr","P","PA","PB","I"]

	validity_syms = filter( kv -> in(kv[1],validity_syms_names), syms)
	symmetries = filter( kv -> in(kv[1],symmetries_names), syms)

	for (k,v) in syms
		@assert in(k,[validity_syms_names;symmetries_names]) "Invalid symmetry name $(k)"
	end



    #get a function to figure out which states are vald
	if constraint === identity
	    is_valid_state = get_valid_state_function(L,unitCellSize,validity_syms)
	else 
		is_valid_state_built_in = get_valid_state_function(L,unitCellSize,validity_syms)
		is_valid_state = x -> (constraint(x) && is_valid_state_built_in(x))
	end


    #figure out which symmetry function we want to use
   	make_orbit = make_easy_orbit_function(L,unitCellSize,symmetries)

    # #count how many basis elements we have
    num_basis_elements :: Int = 0

    #initialize our Dicts
	basis = Dict{UInt64,BASISVECTOR}()
    reps = Dict{UInt64,CONJCLASS}()

    #loop over states in this sector
    for s in 0:(2^L)-1
    	x = UInt64(s)
    	if is_valid_state(x) && !haskey(basis,x)
        	
        	Nullable_orbit = make_orbit(x)


        	if !isnull(Nullable_orbit)
  				orbit = get(Nullable_orbit)

        		#println("state: "*string(digits(x,2,4))*"\t orbit:"*string(orbit))
  				#println(x)
  				#println(orbit)
  				(norm_x,conj_class_x,O_x) = (orbit.norm, orbit.representative, orbit.elements)

  				#we have to correct for the possibility that the representative
  				#doesn't have phase factor 1
  				#in principle clever indexing could prevent that from happening
  				pf_rep = O_x[conj_class_x]
				for (gx, pf_x) in O_x
					#if is_valid_state(gx)
					basis[gx] = BASISVECTOR(conj_class_x, zchop(pf_rep/pf_x))
					#end
				end

                #offset by 1 for 1-indexing
                reps[conj_class_x] =  CONJCLASS(num_basis_elements+1,norm_x)

                #increment number of basis elements
                num_basis_elements += 1
            end
	    end
    end

    if !(constraint === identity)
    	symmetries["constrained"] = 1
    end

    symmetries = merge(symmetries,validity_syms)

	return BASIS(L, unitCellSize, length(reps),basis, reps, symmetries)
end


"""
	check_same_basis_debug(basis1,basis2)

Checks if two bases are the same. UNTESTED! Use for debugging only.
"""
function check_same_basis_debug(
	basis1,
	basis2
	)

    for k in collect(keys(basis1.get_conj_class))
        bv1 = basis1.get_conj_class[k]
        bv2 = basis2.get_conj_class[k]
        if (bv1.conj_class != bv2.conj_class) || (abs(bv1.phase_factor - bv2.phase_factor) >= 10e-5)
            println(k,", ", (bv1.conj_class != bv2.conj_class),", ", abs(bv1.phase_factor - bv2.phase_factor) >= 10e-5)
            println(bv1)
            println(bv2)
        end
    end


    println(length(basis1.conj_classes))
    println(length(basis2.conj_classes))

    for k in collect(keys(basis1.conj_classes))
        cc1 = basis1.conj_classes[k]
        cc2 = basis2.conj_classes[k]
        if cc1.norm != cc2.norm
            println("Bad cc:", k)
            println(cc1)
            println(cc2)
        end
    end
end
