######################################################
#This file uses a conceptually easier way to make basis
#which I anticipate should be slower to use but
#faster to code, which is more important right now.
#
# Requires "basis.jl", since I want to be cross-compatible and
# not re-define things.
#
#################################################



#################### Data Structures ############################


struct BasisVector
    #type for linking a basis vector to it's conjugacy class
    conj_class :: UInt64 ##conjugacy class for the element
    phase_factor :: ComplexF64 ##phase factor to translate
end
Base.show(io::IO, bv::BasisVector) = print(io, "BV[",bin(bv.conj_class), ", ", bv.phase_factor, "]")



struct ConjClass
    #type for describing the conjugacy class of a basis element
    index :: Int64 #index of this element in the Hamiltonian matrix
    norm :: Int #norm^2 of the class so we can normalize it properly
end
Base.show(io::IO, cc::ConjClass) = print(io, "cc[ind:", cc.index,", norm:", cc.norm, "]")


struct Basis
    #type for storing a basis with various quantum numbers
    L :: Int64 #number of sites
    get_conj_class :: Dict{UInt64,BasisVector} #gets basis vector for any state
    conj_classes :: Dict{UInt64,ConjClass}    #gets index and norm of a conjugacy class
    q_numbers :: Dict{String,Int} #list of quantum numbers for the state
    #e.g. q_numbers = [("K",1),("Sz",3)]
end


struct ORBIT
	#type for storing the orbit of a single element under a group action
	norm :: Int
	representative :: UInt64
	elements :: Dict{UInt64,ComplexF64}
end
Base.show(io::IO,orb::ORBIT) = print(io,
"ORBIT[",
"\n\tnorm: ", orb.norm,
"\n\trep: ", orb.representative,
"\n\telements: \n",
[ "$(bin(k)) => $(zchop(v))" for (k,v) in orb.elements],
"]"
	)

##################### Symmetry operations  #########################

function measure_N_up(s::UInt64, L::Int)
    #s - state
    #L - Length of the spin chain, must be even
    #returns - the Sz sector of s, i.e. the number of spins up
    a = Int(sum([(s & (1 << pl)) >> pl for pl = 0:L-1]))
    return a
end 

if testing
    #TEST measure_SZ
    println("Testing N_up")
    f = convert(UInt64, 0)
    a = convert(UInt64, 1)
    b = convert(UInt64, 3)
    c = convert(UInt64, 7)
    d = convert(UInt64, 15)
    display_states([f,a,b,c,d])
    println(map(s -> measure_N_up(s,4),[f,a,b,c,d]))
end


function measure_N_up_A(s :: UInt64, L :: Int)
    return Int(sum([(s & (1 << pl)) >> pl for pl = 0:2:L-1]))
end


if testing
    #TEST measure_SZ
    println("Testing N_up_A")
    f = convert(UInt64, 0)
    a = convert(UInt64, 1)
    b = convert(UInt64, 3)
    c = convert(UInt64, 7)
    d = convert(UInt64, 15)
    display_states([f,a,b,c,d])
    println(map(s -> measure_N_up_A(s,4),[f,a,b,c,d]))
end



function measure_N_up_B(s :: UInt64, L :: Int)
    return Int(sum([(s & (1 << pl)) >> pl for pl = 1:2:L-1]))
end



################ Functions to apply symmetry ###############

"""Returns a function which computes translated versions of states
and gives their phase factors.
"""
function make_translation_function(L :: Int, a :: Int,  K :: Int)
	#L - size of the spin chain
	#a - size of the unit cell
	#K - translation symmetry sector


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

if testing5
	println("Testing 'make_translation_function'. ")
	tr_fcn = make_translation_function(10,2,1)
	println(typeof(tr_fcn))
	x = UInt64(3)
	pf = Complex(1.0)
	println(bin(x),", ", pf)
	(gx,pf2) = tr_fcn(x,pf)
	println(bin(gx),", ", pf2)
end	



"""Returns a function which computes a version of a state flipped on the B sublattice. Needs unitCellSize = 2, of course.
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

# if testing5
# 	println("\nTesting 'make_Z2B_function'. ")
# 	Z2B_fcn = make_Z2B_function(10,1)
# 	Z2B_fcn_2 = make_Z2B_function(10,-1)
# 	x = UInt64(3)
# 	pf = Complex(1.0)
# 	println(bin(x),", ", pf)
# 	(gx,pf2) = Z2B_fcn(x,pf)
# 	println(bin(gx),", ", pf2)
# 	(gx,pf2) = Z2B_fcn_2(x,pf)
# 	println(bin(gx),", ", pf2)
# end




"""Returns a function which computes a version of a state flipped on the A sublattice. Needs unitCellSize = 2, of course.
"""
function make_Z2A_function(
	L :: Int,
	Z2A :: Int
	)

	flipper :: UInt64 = UInt64(sum([1 << k for k in 0:2:L-2]))

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

# if testing5
# 	println("\nTesting 'make_Z2A_function'. ")
# 	Z2A_fcn = make_Z2A_function(10,1)
# 	Z2A_fcn_2 = make_Z2A_function(10,-1)
# 	x = UInt64(3)
# 	pf = Complex(1.0)
# 	println(bin(x),", ", pf)
# 	(gx,pf2) = Z2A_fcn(x,pf)
# 	println(bin(gx),", ", pf2)
# 	(gx,pf2) = Z2A_fcn_2(x,pf)
# 	println(bin(gx),", ", pf2)
# end

"""Returns a function which computes a version of a state flipped on all sites.
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

# if testing6
# 	println("\nTesting 'make_Z2A_function'. ")
# 	Z2A_fcn = make_spin_flip_function(10,1)
# 	Z2A_fcn_2 = make_spin_flip_function(10,-1)
# 	x = UInt64(3)
# 	pf = Complex(1.0)
# 	println(string(x,base=2),", ", pf)
# 	(gx,pf2) = Z2A_fcn(x,pf)
# 	println(string(gx,base=2),", ", pf2)
# 	(gx,pf2) = Z2A_fcn_2(x,pf)
# 	println(string(gx,base=2),", ", pf2)
# end


"""Returns a function which computes a version of a state
	inverted around the middle. Works for any unit cell size.
"""
function make_Inversion_function(
	L :: Int,
	Inv :: Int,
	a :: Int
	)
	
	cell = UInt64((2^a) -1)

	#this seems fairly inefficient
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

if testing5
	println("\nTesting 'make_Inversion_function'. ")
	inv_fcn = make_Inversion_function(10,1,2)
	inv_fcn_2 = make_Inversion_function(10,-1,2)
	x = UInt64(137)
	pf = Complex(1.0)
	println(bin(x),", ", pf)
	(gx,pf2) = inv_fcn(x,pf)
	println(bin(gx),", ", pf2)
	(gx,pf2) = inv_fcn_2(x,pf)
	println(bin(gx),", ", pf2)

end





############ Make a "universal" basis ##################


####symmetries are broken into one of two types:
# 1. symmetries that restrict which states are valid,
# 	 such as magnetization sectors
# 2. and "proper" symmetries with associated phase factors
# The code filters the first out separately, since they're much faster


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
Given the symmetry information, returns a function which computes
the orbit of an individual element. Works for valid basis elements only.
"""
function make_easy_orbit_function(
	L :: Int,
	unitCellSize :: Int,
	symmetries :: Dict{String,Int}
	)
	
	#let's do this the stupid way for now

	has_translation :: Bool = false
	symmetries_fcns = Vector()
	G_size = 1

	#translation is a special case
	if haskey(symmetries,"K")
		sym_fcn = make_translation_function(L,unitCellSize,symmetries["K"])
		sub_gp_size = div(L,unitCellSize)
		push!(symmetries_fcns, (sub_gp_size,sym_fcn))
		has_translation = true
		G_size *= sub_gp_size
	end

	for (name,sector) in symmetries
		if name != "K"
			if name == "Z2A"
				sym_fcn = make_Z2A_function(L,sector)
				sub_gp_size = 2
			elseif name == "Z2B"
				sym_fcn = make_Z2B_function(L,sector)
				sub_gp_size = 2
			elseif name == "Inv"
				sym_fcn = make_Inversion_function(L,sector,unitCellSize)
				sub_gp_size = 2
			elseif name == "Z2"
				sym_fcn = make_spin_flip_function(L,sector)
				sub_gp_size = 2
			end
			G_size *= sub_gp_size
			push!(symmetries_fcns, (sub_gp_size, sym_fcn))
		end
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

if testing5
	println("\nTesting 'make_easy_orbit_function'. ")

	println("\nTesting translation...")
	symmetries = Dict("K" => 1)
	orb_fcn = make_easy_orbit_function(8,2,symmetries)
	
	x = UInt64(1)
	orb = orb_fcn(x)

	println(orb)
	if ~isnull(orb)
		elem = get(orb).elements
		for (k,v) in elem
			println(bin(k), ", ", v)
		end
	end


	println("\nTesting Z2B...")
	symmetries = Dict("Z2B" => 1)
	orb_fcn = make_easy_orbit_function(8,2,symmetries)
	
	x = UInt64(134)
	orb = orb_fcn(x)


	println(orb)

	if ~isnull(orb)
		println("cc_x:", bin(get(orb).representative))
		elem = get(orb).elements
		for (k,v) in elem
			println(bin(k), ", ", v)
		end
	end

	println("\nTesting Z2B + translation...")
	
	symmetries = Dict("Z2B" => -1, "K" => 0)
	orb_fcn = make_easy_orbit_function(8,2,symmetries)
	

	x = UInt64(5)
	orb = orb_fcn(x)

	println(orb)
	if ~isnull(orb)
		elem = get(orb).elements
		for (k,v) in elem
			println(bin(k), ", ", v)
		end
	end



	println("\nTesting Z2B + Z2A +  translation...")
	
	symmetries = Dict("Z2B" => -1, "Z2A" => -1, "K" => 1)
	orb_fcn = make_easy_orbit_function(12,2,symmetries)
	

	x = UInt64(5)
	orb = orb_fcn(x)

	println(orb)
	if ~isnull(orb)
		elem = get(orb).elements
		for (k,v) in elem
			println(bin(k), ", ", v)
		end
	end
end


"""
Produces a basis from its symmetry information.
"""
function make_basis(
	L :: Int;  						#number of lattice sites
	unitCellSize :: Int = 1, 			#self-expanatory
	# validity_syms :: Dict{String,Int} =  Dict{String,Int}(),#symmetries that restrict which states are valid
	# #name => sector, e.g. translation sector 3 is "K" => 3
	syms :: Dict{String,Int} =  Dict{String,Int}(),
	constraint :: Function = identity
	)

	#if we're not using symmetries, then we don't have to compute anything
	#but let's put in placeholders anyway
	if (length(syms) == 0 && constraint === identity)
		return 	Basis(L, 
			Dict{UInt64,BasisVector}(),
			Dict{UInt64,ConjClass}(),
			Dict{String,Int}())
	end



	#separte out the symemtries that simply restrict the Hilbert space
	#hardcoded names for now, as only a few are implemented
	validity_syms_names = ["Sz","SzA","SzB"]
	symmetries_names = ["K","Z2","Z2A","Z2B","Inv"]

	validity_syms = filter( kv -> in(kv[1],validity_syms_names), syms)
	symmetries = filter( kv -> in(kv[1],symmetries_names), syms)

    #get a function to figure out which states are vald
	if constraint === identity
	    is_valid_state = get_valid_state_function(L,unitCellSize,validity_syms)
	else 
		is_valid_state_built_in = get_valid_state_function(L,unitCellSize,validity_syms)
		is_valid_state = x -> (constraint(x,L) && is_valid_state_built_in(x))
	end


    #figure out which symmetry function we want to use
   	make_orbit = make_easy_orbit_function(L,unitCellSize,symmetries)

    # #count how many basis elements we have
    num_basis_elements :: Int = 0

    #initialize our Dicts
	basis = Dict{UInt64,BasisVector}()
    reps = Dict{UInt64,ConjClass}()

    #loop over states in this sector
    for s in 0:(2^L)-1
    	x = UInt64(s)
    	if is_valid_state(x) && !haskey(basis,x)
        	
        	Nullable_orbit = make_orbit(x)


        	if !isnull(Nullable_orbit)
  				orbit = get(Nullable_orbit)
  				#println(x)
  				#println(orbit)
  				(norm_x,conj_class_x,O_x) = (orbit.norm, orbit.representative, orbit.elements)

  				#we have to correct for the possibility that the representative
  				#doesn't have phase factor 1
  				#in principle clever indexing could prevent that from happening
  				pf_rep = O_x[conj_class_x]
				for (gx, pf_x) in O_x
					#if is_valid_state(gx)
					basis[gx] = BasisVector(conj_class_x, zchop(pf_rep/pf_x))
					#end
				end

                #offset by 1 for 1-indexing
                reps[conj_class_x] =  ConjClass(num_basis_elements+1,norm_x)

                #increment number of basis elements
                num_basis_elements += 1
            end
	    end
    end

    if !(constraint === identity)
    	symmetries["constrained"] = 1
    end

	return Basis(L,basis, reps, symmetries)
end


if testing5
    println("Testing easy basis function")

    L = 4

    SzTrZ2basis = make_universal_basis(8,2,Dict("SzA" => 2), Dict("K" => 0))
    println(collect(Set(values(SzTrZ2basis.get_conj_class))))
    println("Basis size:", length(SzTrZ2basis.conj_classes))
    println(SzTrZ2basis.conj_classes)
    display_states(collect(keys(SzTrZ2basis.conj_classes)))
    println(SzTrZ2basis.q_numbers)


    Delta = 0.8
    alpha = 0.1
    g_tau  = 0.3
    u_tau = 0.1
    B_scale = 1.0

    abstract_hamiltonian = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)

    #map(println,make_XXZ_operators(6,0.8))
    H = make_Hamiltonian(L,SzTrZ2basis,abstract_hamiltonian)
    evs = eigfact(Matrix(H))

    println(evs[:values])
end


# ##This is a separate method for now, but it's basically the same thing, so it should be integrated eventually.
# """
# Produces a basis from its symmetry information and a function that constrains the Hilbert space.
# #Arguments
# * 'L :: Int': the number of sites
# * 'unitCellSize :: Int': the number of sites per unit cell
# * 'syms :: Dict{String,Int}': a dictionary of symmetries with their sector names, e.g. the translation sector 3 is "K" => 3
# * 'constraint :: Function': a function (state, L) -> bool to determine if a given state satisfies a constraint
# """
# function make_basis(
# 	L :: Int;  						#number of lattice sites
# 	unitCellSize :: Int = 1, 			#self-expanatory
# 	# validity_syms :: Dict{String,Int} =  Dict{String,Int}(),#symmetries that restrict which states are valid
# 	constraint :: Function = 
# 	syms :: Dict{String,Int} =  Dict{String,Int}(),
# 	)
	
# 	#if we're not using symmetries, then we don't have to compute anything
# 	#but let's put in placeholders anyway
# 	if (length(syms) == 0)
# 		return 	Basis(L, 
# 			Dict{UInt64,BasisVector}(),
# 			Dict{UInt64,ConjClass}(),
# 			Dict{String,Int}())
# 	end

# 	#separte out the symemtries that simply restrict the Hilbert space
# 	#hardcoded names for now, as only a few are implemented
# 	validity_syms_names = ["Sz","SzA","SzB"]
# 	symmetries_names = ["K","Z2A","Z2B","Inv"]

# 	validity_syms = filter( kv -> in(kv[1],validity_syms_names), syms)
# 	symmetries = filter( kv -> in(kv[1],symmetries_names), syms)

#     #get a function to figure out which states are vald
#     is_valid_state_built_in = get_valid_state_function(L,unitCellSize,validity_syms)

#     is_valid_state = x -> (constraint(x,L) && is_valid_state_built_in(x))

#     #figure out which symmetry function we want to use
#    	make_orbit = make_easy_orbit_function(L,unitCellSize,symmetries)

#     # #count how many basis elements we have
#     num_basis_elements :: Int = 0

#     #initialize our Dicts
# 	basis = Dict{UInt64,BasisVector}()
#     reps = Dict{UInt64,ConjClass}()

#     #loop over states in this sector
#     for s in 0:(2^L)-1
#     	x = UInt64(s)
#     	if is_valid_state(x) && !haskey(basis,x)
        	
#         	Nullable_orbit = make_orbit(x)


#         	if !isnull(Nullable_orbit)
#   				orbit = get(Nullable_orbit)
#   				#println(x)
#   				#println(orbit)
#   				(norm_x,conj_class_x,O_x) = (orbit.norm, orbit.representative, orbit.elements)

#   				#we have to correct for the possibility that the representative
#   				#doesn't have phase factor 1
#   				#in principle clever indexing could prevent that from happening
#   				pf_rep = O_x[conj_class_x]
# 				for (gx, pf_x) in O_x
# 					#if is_valid_state(gx)
# 					basis[gx] = BasisVector(conj_class_x, zchop(pf_rep/pf_x))
# 					#end
# 				end

#                 #offset by 1 for 1-indexing
#                 reps[conj_class_x] =  ConjClass(num_basis_elements+1,norm_x)

#                 #increment number of basis elements
#                 num_basis_elements += 1
#             end
# 	    end
#     end
# 	return Basis(L,basis, reps, symmetries)
# end