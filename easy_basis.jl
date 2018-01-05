######################################################
#This file uses a conceptually easier way to make basis
#which I anticipate should be slower to use but
#faster to code, which is more important right now.
#
# Requires "basis.jl", since I want to be cross-compatible and
# not re-define things.
#
#################################################


"""Returns a function which computes translated versions of states
and gives their phase factors.
"""

function make_translation_function(L :: Int, a :: Int,  K :: Int)
	#L - size of the spin chain
	#a - size of the unit cell
	#K - translation symmetry sector

	G_k_size = div(L,a)
	omega :: ComplexF64 = exp(- (2 * pi * im * K)/G_k_size)
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

if testing5
	println("\nTesting 'make_Z2B_function'. ")
	Z2B_fcn = make_Z2B_function(10,1)
	Z2B_fcn_2 = make_Z2B_function(10,-1)
	x = UInt64(3)
	pf = Complex(1.0)
	println(bin(x),", ", pf)
	(gx,pf2) = Z2B_fcn(x,pf)
	println(bin(gx),", ", pf2)
	(gx,pf2) = Z2B_fcn_2(x,pf)
	println(bin(gx),", ", pf2)
end




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

if testing5
	println("\nTesting 'make_Z2A_function'. ")
	Z2A_fcn = make_Z2A_function(10,1)
	Z2A_fcn_2 = make_Z2A_function(10,-1)
	x = UInt64(3)
	pf = Complex(1.0)
	println(bin(x),", ", pf)
	(gx,pf2) = Z2A_fcn(x,pf)
	println(bin(gx),", ", pf2)
	(gx,pf2) = Z2A_fcn_2(x,pf)
	println(bin(gx),", ", pf2)
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
function make_easy_basis(
	L :: Int, 						#number of lattice sites
	unitCellSize :: Int, 			#self-expanatory
	validity_syms :: Dict{String,Int},#symmetries that restrict which states are valid
	#name => sector, e.g. translation sector 3 is "K" => 3
	symmetries :: Dict{String,Int}, #all other symmetires
	)
	

    #get a function to figure out which states are vald
    is_valid_state = get_valid_state_function(L,unitCellSize,validity_symmetries)

    #figure out which symmetry function we want to use
   	make_orbit = make_easy_orbit_function(L,unitCellSize,symmetries)

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
