###########################################################
#This file contains code for producing bases with essentially
#arbitrary symmetries. The code is (or should be) abstract enough to be modular
#and work for pretty much any symmetry.

#Must be used instead of the old case-by-case code

###########################################################




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



############################ Full constructors #################################
#we want a separate constructor for the full basis, since it's so much simpler
#returns a basis object, so it's cross-compatible


function make_full_basis(L :: Int)
    #L - length
    #returns - the full basis of size 2^L, for testing
    basis = Dict{UInt64,BasisVector}([UInt64(s) => BasisVector(UInt64(s),1.0) for s in 0:(2^L)-1])
    #x - > ([x'], phase factor e^i theta(x,x'))
    reps = Dict{UInt64,ConjClass}([UInt64(s) => ConjClass(s+1,1) for s in 0:(2^L)-1])
    #[x] -> (index of [x], Norm([x])^2)
    #offset index by 1 for julia being special
    return Basis(basis, reps, Dict())
end

if testing2
    println("Testing full basis")
    basisFull = make_full_basis(4)
    println(collect(values(basisFull.get_conj_class)))
    display_states(collect(keys(basisFull.conj_classes)))
    println(collect(values(basisFull.conj_classes)))
end 



##################### Symmetry operations  #########################


function apply_T(s::UInt64, L::Int)
    #s - state as a binary integer
    #L - length of the spin chain
    #returns - an Array of UInt64 of translated states
    c1 = UInt64(2^L-1)
    return [((s << g) & c1) | (((s<<g) & ~c1) >> L) for g = 0:(L-1)]
end

if testing
    println("Testing applying translation")
    a = convert(UInt64,1)
    display_states(apply_T(a,4))
end



function apply_T_k(s :: UInt64, a :: Int, L :: Int)
    #s - state as a binary integer
    #L - length of the spin chain
    #a - number of sites in the unit cell
    #returns - an Array of UInt64 of translated states
    c1 = UInt64(2^L-1)
    return [((s << g) & c1) | (((s<<g) & ~c1) >> L) for g = 0: a: (L-1)]
end


if testing
    println("Testing applying translation by a")
    x = convert(UInt64,1)
    display_states(apply_T_k(x,1,6))
    display_states(apply_T_k(x,2,6))
end




function flip_spins(states :: Array{UInt64}, L::Int)
    c = UInt64(sum([1 << k for k in 0:L-1]))
    return map(x -> xor(x,c), states)
end

if testing
    println("Testing spin flipping...")
    a = [UInt64(x) for x in 0:7]
    display_states(a)
    display_states(flip_spins(a,4))
end



function flip_spins_B(states :: Array{UInt64}, L :: Int)
    c = UInt64(sum([1 << (k+1) for k in 0:2:L-2]))
    return map(x -> xor(x,c), states)
end

if testing
    println("Testing spin flipping for B only...")
    a = [UInt64(x) for x in 0:7]
    display_states(a)
    display_states(flip_spins_B(a,6))
end




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

############# Make orbits for the various symmetries ############

###There are two levels of wrappers here right now
#
#The orbit maker calls functions which generate symmetries for individual orbits



function apply_Z2B(L :: Int, Z2 :: Int, old_Orb :: ORBIT)
	cc_x = old_Orb.representative
	norm_x = old_Orb.norm
	O_x = Dict{UInt64,ComplexF64}(old_Orb.elements)
	#bitwise to flip the B spins
	flipper :: UInt64 = UInt64(sum([1 << (k+1) for k in 0:2:L-2]))
	for (gx, pf) in old_Orb.elements

		#compute flipped states
		hgx = xor(gx,flipper)
		
		#this is literally the same code. There should be a better way to do this
		if hgx == old_Orb.representative
			norm_x += pf*Z2*old_Orb.norm
		else
			if !haskey(O_x,hgx)
				O_x[hgx] = pf * Z2
				if hgx < cc_x
					cc_x = hgx
				end
			end
		end
	end
	if abs(norm_x) >= 10e-14
		stab :: Int = length(O_x) == length(old_Orb.elements) ? 2 : 1
		return Nullable(ORBIT(old_Orb.norm * stab , cc_x, O_x))
	else
		return Nullable{ORBIT}()
	end 
end




function apply_K(L :: Int, a :: Int,  K :: Int, x :: UInt64)
	G_k_size = div(L,a)
	w :: ComplexF64 = exp(- (2 * pi * im * K)/G_k_size)
	phase_factor :: ComplexF64 = 1.0
	norm_x :: ComplexF64 = 1.0
	representative :: UInt64 = x

	O_x = Dict{UInt64,ComplexF64}(x => phase_factor)
	sizehint!(O_x,G_k_size)

	c1 = UInt64(2^L-1)
	for g = a : a : (L-1)
		#compute the translated state with fast bitshifts

		shift = x << g
		gx :: UInt64 = (shift & c1) | ((shift & ~c1) >> L)
		phase_factor *= w

		if gx == x
			norm_x += phase_factor
		else 
			if !haskey(O_x,gx) 
				O_x[gx] = phase_factor
				if gx < representative
					representative = gx
				end
			end
		end
	end

	if abs(norm_x) >= 10e-14
		return Nullable(ORBIT(div(G_k_size,length(O_x)), representative, O_x))
	else
		return Nullable{ORBIT}()
	end 
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

###########WARNING!!!!!!!!#######################
#the logic for this function is not final
#and doesn't check that your symmetries are cross-compatible
function get_orbit_function(L :: Int, unitCellSize :: Int, symmetries :: Dict{String,Int})
    #L - length
    #unitCellSize - size of the unit cell
    #symmetries - Array of the symmetries, e.g [("Sz", 2),("K",-2)]
    #returns - a function that creates orbits from valid basis elements

    #this is a placeholder for now
    if haskey(symmetries,"K")
    	if haskey(symmetries,"Z2B")
    		a = unitCellSize
    		Z2B = symmetries["Z2B"]
    		K = symmetries["K"]
    		return function (x)
	    		k_orbit = apply_K(L,a,K,x)
				orbit = !isnull(k_orbit) ? apply_Z2B(L, Z2B, get(k_orbit)) : k_orbit
				return orbit
			end
    	else
    		a = unitCellSize
    		K = symmetries["K"]
    		return function (x)
    			return apply_K(L,a,K,x)
    		end
    	end
    else
    	if haskey(symmetries,"Z2B")
	    	Z2B = symmetries["Z2B"]
	    	return function (x)
	    		return apply_Z2B(L,Z2B,ORBIT(1.0,x,Dict(x => 1.0)))
	    	end
    	else
    		#we have no symmetries, so give the trivial orbit
		    #this is inefficient and I should get a better way to do this
		    return function (x)
		    	return Nullable(ORBIT(1.0,x,Dict(x => 1.0)))
		    end
		end
    end
end



function make_universal_basis(L :: Int, unitCellSize :: Int, validity_symmetries :: Dict{String, Int}, symmetries :: Dict{String,Int})
    #L - length
    #unitCellSize - size of the unit cell
    #symmetries - Array of the symmetries, e.g [("Sz", 2),("K",-2)]

    #get a function to figure out which states are vald
    is_valid_state = get_valid_state_function(L,unitCellSize,validity_symmetries)

    #figure out which symmetry function we want to use
   	make_orbit = get_orbit_function(L,unitCellSize,symmetries)

    # #count how many basis elements we have
    num_basis_elements :: Int = 0

    #initialize our Dicts
	basis = Dict{UInt64,BasisVector}()
    reps = Dict{UInt64,ConjClass}()

    #loop over states in this sector
    for x in [UInt64(s) for s in 0:(2^L)-1]
    	if is_valid_state(x) && !haskey(basis,x)
        	
        	Nullable_orbit = make_orbit(x)

        	if !isnull(Nullable_orbit)
  				orbit = get(Nullable_orbit)
  				(norm_x,conj_class_x,O_x) = (orbit.norm, orbit.representative, orbit.elements)

				for (gx, pf_x) in O_x
					basis[gx] = BasisVector(conj_class_x, numchop(1/pf_x))
				end

                #offset by 1 for 1-indexing
                reps[conj_class_x] =  ConjClass(num_basis_elements+1,norm_x)

                #increment number of basis elements
                num_basis_elements += 1
            end
	    end
    end
	return Basis(basis, reps, symmetries)
end


if testing3
    println("Testing a universal basis function")

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




##########old code for testing


function make_Sz_Tr_flip_basis(L :: Int, N_up_A :: Int, K :: Int, Z2B_flip :: Int)
    #L - length
    #N_up - number of spin ups
    #k - translation sector
    #Z2_flip - +1/-1 for flipping the state
    # if N_up_A == div(L,4)
    #     error("This requires spin inversion for A also")
    # end

    #size of the unit cell
    a = 2
    #size of the group
    G_size = L
    #size of the translation part of the group
    G_k_size = div(L,2)

    SzSector = filter(s -> measure_N_up_A(s,L) == N_up_A, [UInt64(s) for s in 0:(2^L)-1])

    #really this should be some sorted data structure, but
    #not sure how Julia stores things internally so
    #a dict for now
    elements = Dict{UInt64,UInt64}()
    #count how many basis elements we have
    num_basis_elements = 0

    #initialize our arrays
    basis_array = Array{Tuple{UInt64,BasisVector},1}(uninitialized,0)
    reps_array = Array{Tuple{UInt64,ConjClass},1}(uninitialized,0)
    

    #loop over states in this sector
    for x in SzSector
        #if x is not present in elements, do...
        get(elements,x) do
            translated_states = apply_T_k(x,a,L)
            #display_states(translated_states)
            k_orbit_size = length(Set(translated_states))

            # #only proceed if we have a non-zero norm
            if mod(K * k_orbit_size,G_k_size) == 0
                #now compute the flipped states too
                orbit_x = [translated_states; flip_spins_B(translated_states,L)]
                orbit_size_x = length(Set(orbit_x))

                #compute our phase factors
                phase_factors_k = [exp(- (2 * pi * im * K * r)/G_k_size) for r in 0 : G_k_size-1]
                phase_factors = [phase_factors_k; Z2B_flip*phase_factors_k]

                stabilizer = filter(r -> orbit_x[r] == x, 1:L)

                norm_x = numchop(sum([phase_factors[r] for r in stabilizer]))


                if abs(norm_x) > 10e-14

                    if abs(norm_x - (G_size/orbit_size_x)) > 10e-10
                        error("error: Bad norm size", norm_x, " orbit_size: ", orbit_size_x)
                    end
                    #find the new representative state
                    conj_class_x = minimum(orbit_x)

                    elem = Dict([y => conj_class_x for y in orbit_x])
                    merge!(elements,elem)

                    basis_array_x = collect(Set([(orbit_x[r],
                                     BasisVector(conj_class_x,numchop(1/phase_factors[r])))
                                     for r in 1:G_size]))
                    append!(basis_array,basis_array_x)

                    if imag(norm_x) >= 10e-15
                        error("error: Large imaginary part of norm!")
                    end
                    #offset by 1 for 1-indexing
                    reps_array_x = (conj_class_x, ConjClass(num_basis_elements+1,real(norm_x)))
                    push!(reps_array,reps_array_x)

                    #increment number of basis elements
                    num_basis_elements += 1
                end
            end
        end
    end

    basis = Dict{UInt64,BasisVector}(basis_array)
    reps = Dict{UInt64,ConjClass}(reps_array)
    sym  = Dict("N_up_A"=>N_up_A,"K"=>K,"Z_2^B (flip)"=>Z2B_flip)
    return Basis(basis, reps, sym)
end


function make_Sz_Tr_basis(L :: Int, N_up_A :: Int, K :: Int)
    #L - length
    #N_up - number of spin ups
    #k - translation sector

    #non-trivial unit-cell size
    a = 2 

    SzSector = filter(s -> measure_N_up_A(s,L) == N_up_A, [UInt64(s) for s in 0:(2^L)-1])

    # println("Sz length: ", length(SzSector))


    #really this should be some sorted data structure, but
    #not sure how Julia stores things internally so
    #a dict for now
    elements = Dict{UInt64,UInt64}()
    #count how many basis elements we have
    num_basis_elements = 0

    #initialize our arrays
    basis_array = Array{Tuple{UInt64,BasisVector},1}(uninitialized,0)
    reps_array = Array{Tuple{UInt64,ConjClass},1}(uninitialized,0)
    
    #size of the group
    G_size = div(L,2)
    #size of the translation part of the group
    G_k_size = G_size

    #loop over states in this sector
    for x in SzSector
        #if x is not present in elements, do...
        get(elements,x) do
            translated_states = apply_T_k(x,a,L)
            #display_states(translated_states)
            k_orbit_size = length(Set(translated_states)) 

            #only proceed if we have a non-zero norm
            if mod(K * k_orbit_size,G_k_size) == 0
                conj_class_x = minimum(translated_states)

                elem = Dict([y => conj_class_x for y in translated_states])
                merge!(elements,elem)

                #compute our phase factors
                phase_factors = [exp(- (2 * pi * im * K * r)/G_k_size) for r in 0 : k_orbit_size-1]

                k_norm = G_size/k_orbit_size

                basis_array_x = [(translated_states[r],
                                 BasisVector(conj_class_x,numchop(1/phase_factors[r])))
                                 for r in 1:k_orbit_size]
                append!(basis_array,basis_array_x)

                if imag(k_norm) >= 10e-15
                    error("WARNING: Large imaginary part of norm!")
                end
                #offset by 1 for 1-indexing
                reps_array_x = (conj_class_x, ConjClass(num_basis_elements+1,real(k_norm)))
                push!(reps_array,reps_array_x)

                #increment number of basis elements
                num_basis_elements += 1
            end
        end
    end

    basis = Dict{UInt64,BasisVector}(basis_array)
    reps = Dict{UInt64,ConjClass}(reps_array)
    sym  = Dict("N_up_A"=>N_up_A,"K"=>K)
    return Basis(basis, reps, sym)
end