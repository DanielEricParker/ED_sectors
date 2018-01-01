######################################################################
#produces basis for a given symmetry sector
#right now this is all case-by-case

#I think case-by-case is actually necessary to make everything the most efficient
#but it's not very satisfying
#a more abstract approach would be nice eventually

#####################################################################


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
    q_numbers :: Array{Tuple{String,Int}} #list of quantum numbers for the state
    #e.g. q_numbers = [("K",1),("Sz",3)]
end

#more general version for the "univeral" constructor
struct Basis2
    #type for storing a basis with various quantum numbers
    get_conj_class :: Dict{UInt64,BasisVector} #gets basis vector for any state
    conj_classes :: Dict{UInt64,ConjClass}    #gets index and norm of a conjugacy class
    q_numbers :: Dict{String,Int} #list of quantum numbers for the state
    #e.g. q_numbers = [("K",1),("Sz",3)]
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






############################ Basis constructors #################################


function make_full_basis(L :: Int)
    #L - length
    #returns - the full basis of size 2^L, for testing
    basis = Dict{UInt64,BasisVector}([UInt64(s) => BasisVector(UInt64(s),1.0) for s in 0:(2^L)-1])
    #x - > ([x'], phase factor e^i theta(x,x'))
    reps = Dict{UInt64,ConjClass}([UInt64(s) => ConjClass(s+1,1) for s in 0:(2^L)-1])
    #[x] -> (index of [x], Norm([x])^2)
    #offset index by 1 for julia being special
    return Basis(basis, reps,[])
end

if testing2
    println("Testing full basis")
    basisFull = make_full_basis(4)
    println(collect(values(basisFull.get_conj_class)))
    display_states(collect(keys(basisFull.conj_classes)))
    println(collect(values(basisFull.conj_classes)))
end 



################## Sz sectors basis  ###############3

function make_Sz_basis(L :: Int, N_up_A :: Int)
    #L - length
    #N_up - number of spin up's in the sector we care about
    #returns - the full basis of size 2^L, for testing

    SzSector = filter(s -> measure_N_up_A(s,L) == N_up_A, [UInt64(s) for s in 0:(2^L)-1])
    display_states(SzSector)
    #display_states(SzSector)
    basis = Dict{UInt64,BasisVector}([UInt64(s) => BasisVector(s,1.0) for s in SzSector])
    #x - > ([x'], phase factor e^i theta(x,x'))
    reps = Dict{UInt64,ConjClass}([UInt64(s) => ConjClass(index,1) for (index, s) in enumerate(SzSector)])
    #[x] -> (index of [x], Norm([x])^2)
    #offset index by 1 for julia being special
    return Basis(basis, reps, [("N_up_A", N_up_A)])
end

if testing
    println("Testing making a basis with Sz sector for XXZ*" )
    println(make_Sz_basis(4,3))

    L = 8
    basisSz = make_Sz_basis(L,2)

    Delta = 0.8
    alpha = 0.1
    g_tau  = 0.3
    u_tau = 0.1
    B_scale = 1.0

    abstract_hamiltonian = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)
    #map(println,make_XXZ_operators(6,0.8))
    H = make_Hamiltonian(L,basisSz,abstract_hamiltonian)
    evs = eigfact(Matrix(H))

    println(evs[:values][1:40])

    #this seems to work against quspin

end


############### Sz and k sectors basis ################

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

    return Basis(basis, reps,[("N_up_A",N_up_A),("K",K)])
end

if testing3
    println("Testing a basis with both Sz and Translation")

    L = 8

    SzTrbasis = make_Sz_Tr_basis(8,2,0)
    println("Basis size:", length(SzTrbasis.conj_classes))
    println(collect(Set(values(SzTrbasis.get_conj_class))))


    # Delta = 0.8
    # alpha = 0.1
    # g_tau  = 0.3
    # u_tau = 0.1
    # B_scale = 1.0

    # abstract_hamiltonian = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)
    # #map(println,make_XXZ_operators(6,0.8))
    # H = make_Hamiltonian(L,SzTrbasis,abstract_hamiltonian)
    # evs = eigfact(Matrix(H))

    # println(evs[:values])

    #println(SzTrbasis.get_conj_class[UInt(153)])

    #this seems to work against quspin
end




############ Make a basis with Sz, Translation, and Z_2 (flip) sectors

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

    return Basis(basis, reps,[("N_up_A",N_up_A),("K",K),("Z_2^B (flip)",Z2B_flip)])
end


if testing
    println("Testing a basis with, Sz, Translation and Z2 (flip)")

    L = 8

    SzTrZ2basis = make_Sz_Tr_flip_basis(8,2,0,-1)
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




############ Make a "universal" basis for any implemented symmetry


function get_orbit_function(L :: Int, unitCellSize :: Int, symmetries :: Dict{String,Int})
    #L - length
    #unitCellSize - size of the unit cell
    #symmetries - Array of the symmetries, e.g [("Sz", 2),("K",-2)]
    #returns - a function that creates orbits from valid basis elements

    #this is a placeholder for now

    return x -> make_orbit_K(L, unitCellSize, symmetries["K"],x)
end


struct ORBIT
	norm :: Float64
	representative :: UInt64
	elements :: Dict{UInt64,ComplexF64}
end

function make_orbit_K(L :: Int, a :: Int, K :: Int, x :: UInt64)

	G_k_size = div(L,a)
	w = exp(- (2 * pi * im * K)/G_k_size)
	phase_factor :: ComplexF64 = 1.0
	norm_x :: ComplexF64 = 1.0
	representative = x

	O_x = Dict{UInt64,ComplexF64}(x => phase_factor)
	sizehint!(O_x,G_k_size)

	c1 = UInt64(2^L-1)
	for g = a : a : (L-1)
		shift = x << g
		gx = (shift & c1) | ((shift & ~c1) >> L)
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
		return Nullable(ORBIT(real(numchop(norm_x)), representative, O_x))
	else
		return Nullable{ORBIT}()
	end 
end


if testing3
	println("Testing making translation orbit")
	orb = get(make_orbit_K(4,1,1,UInt64(1)))
	println(orb)
	k = collect(keys(orb.elements))
	display_states(k)
end

function make_universal_basis(L :: Int, unitCellSize :: Int, symmetries :: Dict{String,Int})
    #L - length
    #unitCellSize - size of the unit cell
    #symmetries - Array of the symmetries, e.g [("Sz", 2),("K",-2)]

    #list all our states. For *really* huge sizes this is impractical, but it's easier here for now
    states = [UInt64(s) for s in 0:(2^L)-1]

    #if we have an SzA sector restriction, use it first
    if haskey(symmetries,"SzA") 
    	N_up_A = symmetries["SzA"]
    	states = filter(s -> measure_N_up_A(s,L) == N_up_A, states)
    end

    #figure out which symmetry function we want to use
   	make_orbit = get_orbit_function(L,unitCellSize,symmetries)

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
    for x in states
        #if x is not present in elements, do...
        get(elements,x) do
        	
        	orbit = make_orbit(x)

        	if !isnull(orbit)
  				orbit = get(orbit)
  				(norm_x,conj_class_x,O_x) = (orbit.norm, orbit.representative, orbit.elements)

                elem = Dict([y => conj_class_x for y in collect(keys(O_x))])
                merge!(elements,elem)

                # basis_array_x = collect(Set([(orbit_x[r],
                #                  BasisVector(conj_class_x,numchop(1/phase_factors[r])))
                #                  for r in 1:G_size]))
                basis_array_x = [(gx, BasisVector(conj_class_x, numchop(1/pf_x))) for (gx, pf_x) in O_x]
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

    basis = Dict{UInt64,BasisVector}(basis_array)
    reps = Dict{UInt64,ConjClass}(reps_array)

    #return Basis2(basis, reps,symmetries])
	return Basis(basis, reps, [Tuple(kv) for kv in symmetries])
end


if testing3
    println("Testing a universal basis function")

    L = 4

    SzTrZ2basis = make_universal_basis(8,2,Dict("SzA" => 2,"K" => 0))
    println(collect(Set(values(SzTrZ2basis.get_conj_class))))
    println("Basis size:", length(SzTrZ2basis.conj_classes))
    println(SzTrZ2basis.conj_classes)
    display_states(collect(keys(SzTrZ2basis.conj_classes)))
    println(SzTrZ2basis.q_numbers)


    # Delta = 0.8
    # alpha = 0.1
    # g_tau  = 0.3
    # u_tau = 0.1
    # B_scale = 1.0

    # abstract_hamiltonian = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)

    # #map(println,make_XXZ_operators(6,0.8))
    # H = make_Hamiltonian(L,SzTrZ2basis,abstract_hamiltonian)
    # evs = eigfact(Matrix(H))

    # println(evs[:values])
end
