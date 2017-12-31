#!/usr/bin/env julia

##################### are we testing? #######################

testing = false
testing2 = true




#################### Display/pretty printing ##########################

function display_states(s::Array{UInt64})
    #x - state
    #returns - nothing
    #displays
    println(map(bin,s))
end


################### General little useful functions #############


function numchop(x :: Complex{Float64}, epsilon = 10e-15)
    #x - complex number
    #returns - x with small real and imaginary parts removed
    a = real(x)
    b = imag(x)

    w = abs(a) < epsilon ? 0 : a
    z = abs(b) < epsilon ? 0 : b 

    return complex(w,z)
end


#################### Data Structures ############################


struct VEC
    #type for storing factor * basis_state
    factor :: Float64
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
end
Base.show(io::IO, t::TERM) = print(io,  t.prefactor, "*", t.operator)



function term(prefactor :: Float64, ops :: OP ...)
	#alternative constructor for TERM type
	#uses variable arguments for improved readability
	return TERM(prefactor, collect(ops))
end

if testing2
    testTerm = TERM(0.5,[OP("+",3),OP("-",4),OP("-",5)])
    println(testTerm)
    println(typeof(testTerm.operator))
end

if testing2
	println("Testing making term in the Hamiltonian")
	println(term(0.5,OP("+",3),OP("-",4),OP("-",5)))
end

struct HAMILTONIAN
    #type for storing Hamiltonians
    name :: String
    terms :: Array{TERM}
end

function Base.show(io::IO, H::HAMILTONIAN)
	print(io, "Hamiltonian: ",  H.name, " (Length ", length(H.terms), ")")
	for t in H.terms
		print(io,"\n", t)
	end
end


function Base.:+(H :: HAMILTONIAN, t :: TERM)
	#gives a way to add terms to Hamiltonians
	#usage: h += term
	terms = H.terms
	push!(terms,t)
	HAMILTONIAN(H.name,terms)
end

if testing2
    println("Testing making a Hamiltonian")
    t1 = term(0.5,OP("+",3),OP("-",4),OP("-",5))
    t2 = term(0.5,OP("-",3),OP("z",5))
    println(typeof(t1.operator))
    println(typeof(t2.operator))
    testHam = HAMILTONIAN("testing", [t1])
    println(testHam)
    
    testHam += t2
    println(testHam)
end


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


function apply_operators(v :: VEC, prefactor :: Float64, ops :: Array{OP})
    # v vector
    # list of operators to apply to various sites
    # only works for operators whose output is a single basis vector
    for op in ops
        if is_Null_VEC(v)
            return v
        else 
            if op.name == "X"
                v = apply_S_x(v,op.site)
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
    v = VEC(v.factor * prefactor,v.state)
    return v
end

if testing2
    ops = [OP("X",1),OP("X",2)]
    ops4 = [OP("-",5)]
    v = VEC(1.0,15)
    println("Testing Applying Multiple Operators")
    println(v)
    map(println,ops)
    println(apply_operators(v,2.0,ops))
    println(apply_operators(v,3.0,ops4))
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

##################### Make the operators for the XXZ Hamiltonian ##############



function make_XXZ_operators(L :: Int, Delta :: Float64)
    #L -- length of the spin chain
    #Delta -- anisotropy parameter
    opsList = []
    for i in 0:L-1
        push!(opsList, (0.5, [OP("+",i),OP("-",mod(i+1,L))]))
        push!(opsList, (0.5, [OP("-",i),OP("+",mod(i+1,L))]))
        push!(opsList, (Delta*1.0, [OP("Z",i),OP("Z",mod(i+1,L))]))
    end

    return opsList
end

if testing
    println("Testing Making XXZ Operators")
    map(println,make_XXZ_operators(6,0.8))
end




function make_XXZ_operators_new(L :: Int, Delta :: Float64)
    #L -- length of the spin chain
    #Delta -- anisotropy parameter
    H = HAMILTONIAN("XXZ", [])
    for i in 0:L-1
        H += term(0.5, OP("+",i), OP("-",(i+1)%L))
        H += term(0.5, OP("-",i), OP("+",(i+1)%L))
        H += term(Delta*1.0, OP("Z",i), OP("Z",(i+1)%L))
    end

    return H
end

if testing
	println(make_XXZ_operators_new(4,0.7))
end


function make_XXZ_star_operators(L :: Int, Delta :: Float64, alpha :: Float64, g_tau :: Float64, u_tau :: Float64, B_scale :: Float64)

    opsList = []
    for i in 0:2:L-1
        push!(opsList, (0.5, [OP("+",i),OP("-",mod(i+2,L))]))
        push!(opsList, (0.5, [OP("-",i),OP("+",mod(i+2,L))]))
        push!(opsList, (Delta*1.0, [OP("Z",i),OP("Z",mod(i+2,L))]))
        push!(opsList, (-1.0*B_scale, [OP("X",mod(i+1,L))]))
        push!(opsList, (-g_tau*B_scale, [OP("Z",mod(i+1,L)),OP("Z",mod(i+3,L))]))
        push!(opsList, (-u_tau*B_scale, [OP("X",mod(i+1,L)),OP("X",mod(i+3,L))]))
        push!(opsList, (alpha*1.0, [OP("Z",i),OP("X",mod(i+1,L)), OP("Z",mod(i+2,L))]))
    end
    return opsList
end


if testing2
    println("Testing abstract basis for XXZ*")
    L=4
    Delta = 0.8
    alpha = 0.0
    g_tau  = 0.3
    u_tau = 0.1
    B_scale = 1.0
    H_XXZ_1 = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)

end



function make_XXZ_star_operators_new(L :: Int, Delta :: Float64, alpha :: Float64, g_tau :: Float64, u_tau :: Float64, B_scale :: Float64)

	H = HAMILTONIAN("XXZ_star",[])
    for i in 0:2:L-1
		H += term(0.5, 				OP("+", i), 	OP("-", (i+2) % L))
		H += term(0.5, 				OP("-", i), 	OP("+", (i+2) % L))
		H += term(Delta, 			OP("Z", i), 	OP("Z", (i+2) % L))
	 	H += term(-1.0*B_scale, 	OP("X",(i+1) % L))
	 	H += term(-g_tau*B_scale, 	OP("Z",(i+1) % L),OP("Z",(i+3) % L))
	 	H += term(-u_tau*B_scale, 	OP("X",(i+1) % L),OP("X",(i+3) % L))
	 	H += term(alpha, 			OP("Z",i),		OP("X",(i+1) % L),	OP("Z", (i+2)%L))
    end
    return H
end

if testing2
    println("Testing abstract basis for XXZ*")
    L=4
    Delta = 0.8
    alpha = 0.0
    g_tau  = 0.3
    u_tau = 0.1
    B_scale = 1.0
    H_XXZ_2 = make_XXZ_star_operators_new(L,Delta,alpha,g_tau,u_tau,B_scale)
    [println(TERM(H_XXZ_1[i][1],H_XXZ_1[i][2]),"\t", H_XXZ_2.terms[i]) for i in 9:10]
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


################ Hamiltonian constructor for a given basis ##########################


function make_Hamiltonian(L :: Int, basis :: Basis, opsList :: Array{Any})
    #L - length of the spin chain
    #basis - basis Dict
    #reps - represenatives of conjugacy classes Dict
    #opsList - list of operators we want to add to the Hamiltonian


    dim = length(basis.conj_classes)
    H = spzeros(ComplexF64,dim,dim)

    #iterate over conjugacy classes
    for x in keys(basis.conj_classes)
        cc_x = basis.conj_classes[x]

        #loop over terms in the Hamiltonian
        for (prefactor, ops) in opsList
            #find |y> = H |x>
            y = apply_operators(VEC(1.0,x),prefactor,ops)

            if ~is_Null_VEC(y)
                bv_y = get(basis.get_conj_class,y.state,-1)

                #check if the element exists --- it could have zero norm
                if bv_y != -1
                    cc_y = basis.conj_classes[bv_y.conj_class]

                    #check if upper diagonal
                    if cc_y.index >= cc_x.index
                        #compute the matrix element --- see Note
                        h_xy = sqrt(cc_y.norm/cc_x.norm) * bv_y.phase_factor * y.factor
                        H[cc_x.index,cc_y.index] += h_xy

                        #check if diagonal
                        if cc_y.index > cc_x.index
                             H[cc_y.index,cc_x.index] += conj(h_xy)
                        end
                    end
                # else
                #     println("Warning: no target in basis for state ",y.state)
                end
            end
        end
    end
    return H
end 

if testing2
    println("Testing making the Hamiltonian from the full basis")
    abstract_hamiltonian = make_XXZ_operators(4,0.8)
    #map(println,make_XXZ_operators(6,0.8))
    H1 = make_Hamiltonian(4,basisFull,abstract_hamiltonian)
end



function make_Hamiltonian(L :: Int, basis :: Basis, abstract_Ham :: HAMILTONIAN)
    #L - length of the spin chain
    #basis - basis Dict
    #reps - represenatives of conjugacy classes Dict
    #opsList - list of operators we want to add to the Hamiltonian


    dim = length(basis.conj_classes)
    H = spzeros(ComplexF64,dim,dim)

    #iterate over conjugacy classes
    for x in keys(basis.conj_classes)
        cc_x = basis.conj_classes[x]

        #loop over terms in the Hamiltonian
        for term in abstract_hamiltonian.terms
            #find |y> = H |x>
            y = apply_operators(VEC(1.0,x),term)

            if ~is_Null_VEC(y)
                bv_y = get(basis.get_conj_class,y.state,-1)

                #check if the element exists --- it could have zero norm
                if bv_y != -1
                    cc_y = basis.conj_classes[bv_y.conj_class]

                    #check if upper diagonal
                    if cc_y.index >= cc_x.index
                        #compute the matrix element --- see Note
                        h_xy = sqrt(cc_y.norm/cc_x.norm) * bv_y.phase_factor * y.factor
                        H[cc_x.index,cc_y.index] += h_xy

                        #check if diagonal
                        if cc_y.index > cc_x.index
                             H[cc_y.index,cc_x.index] += conj(h_xy)
                        end
                    end
                # else
                #     println("Warning: no target in basis for state ",y.state)
                end
            end
        end
    end
    return H
end 

if testing2
    println("Testing making the Hamiltonian from the full basis")
    abstract_hamiltonian = make_XXZ_operators_new(4,0.8)
   	#println(abstract_hamiltonian)
    H2 = make_Hamiltonian(4,basisFull,abstract_hamiltonian)
    println(H1 == H2)
end



# #################### Symmetry operations  #########################


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



# # ################### test against quspin #########
# L = 8


# basisFull = make_full_basis(L)

# Delta = 0.8
# alpha = 0.1
# g_tau  = 0.3
# u_tau = 0.1
# B_scale = 1.0

# abstract_hamiltonian = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)
# #map(println,make_XXZ_operators(6,0.8))
# H = make_Hamiltonian(L,basisFull,abstract_hamiltonian)
# evs = eigfact(Matrix(H))

# println(evs[:values][1:40])
# #this works against quspin, so I think we can use this for a baseline against other sectors
# #with XXZ*




# ################## Sz sectors basis  ###############3

# function make_Sz_basis(L :: Int, N_up_A :: Int)
#     #L - length
#     #N_up - number of spin up's in the sector we care about
#     #returns - the full basis of size 2^L, for testing

#     SzSector = filter(s -> measure_N_up_A(s,L) == N_up_A, [UInt64(s) for s in 0:(2^L)-1])
#     display_states(SzSector)
#     #display_states(SzSector)
#     basis = Dict{UInt64,BasisVector}([UInt64(s) => BasisVector(s,1.0) for s in SzSector])
#     #x - > ([x'], phase factor e^i theta(x,x'))
#     reps = Dict{UInt64,ConjClass}([UInt64(s) => ConjClass(index,1) for (index, s) in enumerate(SzSector)])
#     #[x] -> (index of [x], Norm([x])^2)
#     #offset index by 1 for julia being special
#     return Basis(basis, reps, [("N_up_A", N_up_A)])
# end

# # if testing
# #     println("Testing making a basis with Sz sector for XXZ*" )
# #     println(make_Sz_basis(4,3))

# #     L = 8
# #     basisSz = make_Sz_basis(L,2)

# #     Delta = 0.8
# #     alpha = 0.1
# #     g_tau  = 0.3
# #     u_tau = 0.1
# #     B_scale = 1.0

# #     abstract_hamiltonian = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)
# #     #map(println,make_XXZ_operators(6,0.8))
# #     H = make_Hamiltonian(L,basisSz,abstract_hamiltonian)
# #     evs = eigfact(Matrix(H))

# #     println(evs[:values][1:40])

# #     #this seems to work against quspin

# # end


# ############### Sz and k sectors basis ################

# function make_Sz_Tr_basis(L :: Int, N_up_A :: Int, K :: Int)
#     #L - length
#     #N_up - number of spin ups
#     #k - translation sector

#     #non-trivial unit-cell size
#     a = 2 

#     SzSector = filter(s -> measure_N_up_A(s,L) == N_up_A, [UInt64(s) for s in 0:(2^L)-1])

#     # println("Sz length: ", length(SzSector))


#     #really this should be some sorted data structure, but
#     #not sure how Julia stores things internally so
#     #a dict for now
#     elements = Dict{UInt64,UInt64}()
#     #count how many basis elements we have
#     num_basis_elements = 0

#     #initialize our arrays
#     basis_array = Array{Tuple{UInt64,BasisVector},1}(uninitialized,0)
#     reps_array = Array{Tuple{UInt64,ConjClass},1}(uninitialized,0)
    
#     #size of the group
#     G_size = div(L,2)
#     #size of the translation part of the group
#     G_k_size = G_size

#     #loop over states in this sector
#     for x in SzSector
#         #if x is not present in elements, do...
#         get(elements,x) do
#             translated_states = apply_T_k(x,a,L)
#             #display_states(translated_states)
#             k_orbit_size = length(Set(translated_states)) 

#             #only proceed if we have a non-zero norm
#             if mod(K * k_orbit_size,G_k_size) == 0
#                 conj_class_x = minimum(translated_states)

#                 elem = Dict([y => conj_class_x for y in translated_states])
#                 merge!(elements,elem)

#                 #compute our phase factors
#                 phase_factors = [exp(- (2 * pi * im * K * r)/G_k_size) for r in 0 : k_orbit_size-1]

#                 k_norm = G_size/k_orbit_size

#                 basis_array_x = [(translated_states[r],
#                                  BasisVector(conj_class_x,numchop(1/phase_factors[r])))
#                                  for r in 1:k_orbit_size]
#                 append!(basis_array,basis_array_x)

#                 if imag(k_norm) >= 10e-15
#                     error("WARNING: Large imaginary part of norm!")
#                 end
#                 #offset by 1 for 1-indexing
#                 reps_array_x = (conj_class_x, ConjClass(num_basis_elements+1,real(k_norm)))
#                 push!(reps_array,reps_array_x)

#                 #increment number of basis elements
#                 num_basis_elements += 1
#             end
#         end
#     end

#     basis = Dict{UInt64,BasisVector}(basis_array)
#     reps = Dict{UInt64,ConjClass}(reps_array)

#     return Basis(basis, reps,[("N_up_A",N_up_A),("K",K)])
# end

# if testing
#     println("Testing a basis with both Sz and Translation")

#     L = 8

#     SzTrbasis = make_Sz_Tr_basis(8,2,0)
#     println("Basis size:", length(SzTrbasis.conj_classes))
#     println(collect(Set(values(SzTrbasis.get_conj_class))))


#     Delta = 0.8
#     alpha = 0.1
#     g_tau  = 0.3
#     u_tau = 0.1
#     B_scale = 1.0

#     abstract_hamiltonian = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)
#     #map(println,make_XXZ_operators(6,0.8))
#     H = make_Hamiltonian(L,SzTrbasis,abstract_hamiltonian)
#     evs = eigfact(Matrix(H))

#     println(evs[:values])

#     println(SzTrbasis.get_conj_class[UInt(153)])

#     #this seems to work against quspin
# end



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


################ All sectors ###################

#generate all eigenvalues in all sectors to compare with full ED


# L = 20

# Delta = 0.8
# alpha = 0.1
# g_tau  = 0.3
# u_tau = 0.1
# B_scale = 1.0
# eigenvalues = []
# for Sz in 0:div(L,2)
#     for K in 0:div(L,2)-1
#         for Z2B in [-1,1]
#             SzTrZ2basis = make_Sz_Tr_flip_basis(L,Sz,K,Z2B)
#             println("Sz: ", Sz, " K: ", K, " Z2B: ", Z2B," Length; ", length(collect(Set(values(SzTrZ2basis.conj_classes)))))
#             H = make_Hamiltonian(L,SzTrZ2basis,abstract_hamiltonian)
#             evs = eigfact(Matrix(H))
#             max_imag_part = maximum(abs(imag(evs[:values])))
#             if(max_imag_part >= 10e-6)
#                 error("Large imaginary part:", max_imag_part)
#             end
#             append!(eigenvalues,evs[:values])
#         end
#     end
# end
# println("Number of eigenvalues: ", length(eigenvalues), ". Should be: ", 2^L)

# eigenvalues2 = sort(real(eigenvalues))
# writedlm("data/sectors_eigenvalues_$(L)_alpha_$(alpha).csv",eigenvalues2,", ")






# L = 20
# Sz = div(L,4)+3
# K = 3
# Z2B = 1
# println("Making basis...")
# SzTrZ2basis = make_Sz_Tr_flip_basis(L,Sz,K,Z2B)
# println("Made basis for Sz: ", Sz, " K: ", K, " Z2B: ", Z2B," Length: ", length(collect(Set(values(SzTrZ2basis.conj_classes)))))
# for alpha in 0.0 : 0.05 : 0.35
#     println("alpha: ", alpha)
#     Delta = -0.8
#     #alpha = 0.1
#     g_tau  = 0.3
#     u_tau = 0.1
#     B_scale = 1.0
#     abstract_hamiltonian = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)

#     println("Making Hamiltonian...")
#     H = make_Hamiltonian(L,SzTrZ2basis,abstract_hamiltonian)
#     println("Diagonalizing Hamiltonian...")
#     evs = eigfact(Matrix(H))
#     max_imag_part = maximum(abs.(imag(evs[:values])))
#     if(max_imag_part >= 10e-6)
#         error("Large imaginary part:", max_imag_part)
#     end

#     evs2 = real(evs[:values])

#     file = @sprintf "data/sectors_eigenvalues_%i_SzA_%i_K_%i_Z2B_%i_delta_%.4f_alpha_%.4f.csv" L Sz K Z2B Delta alpha
#     println(file)
#     writedlm(file,evs2,", ")
# end