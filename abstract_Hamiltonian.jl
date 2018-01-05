#####################################################################
#This file gives tool for constructing fairly arbitrary spin Hamiltonians

#The main tool is an "abstract Hamiltonian" (i.e. human-readable)
#type called Hamiltonian

#and a way to convert HAMILTONIAN + BASIS into an actual matrix

#####################################################################



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



#########################Example Hamiltonians ########################



# function make_XXZ_operators(L :: Int, Delta :: Float64)
#     #L -- length of the spin chain
#     #Delta -- anisotropy parameter
#     opsList = []
#     for i in 0:L-1
#         push!(opsList, (0.5, [OP("+",i),OP("-",mod(i+1,L))]))
#         push!(opsList, (0.5, [OP("-",i),OP("+",mod(i+1,L))]))
#         push!(opsList, (Delta*1.0, [OP("Z",i),OP("Z",mod(i+1,L))]))
#     end

#     return opsList
# end

# if testing
#     println("Testing Making XXZ Operators")
#     map(println,make_XXZ_operators(6,0.8))
# end




# function make_XXZ_operators_new(L :: Int, Delta :: Float64)
#     #L -- length of the spin chain
#     #Delta -- anisotropy parameter
#     H = HAMILTONIAN("XXZ", [])
#     for i in 0:L-1
#         H += term(0.5, OP("+",i), OP("-",(i+1)%L))
#         H += term(0.5, OP("-",i), OP("+",(i+1)%L))
#         H += term(Delta*1.0, OP("Z",i), OP("Z",(i+1)%L))
#     end

#     return H
# end

# if testing
# 	println(make_XXZ_operators_new(4,0.7))
# end


# function make_XXZ_star_operators(L :: Int, Delta :: Float64, alpha :: Float64, g_tau :: Float64, u_tau :: Float64, B_scale :: Float64)

#     opsList = []
#     for i in 0:2:L-1
#         push!(opsList, (0.5, [OP("+",i),OP("-",mod(i+2,L))]))
#         push!(opsList, (0.5, [OP("-",i),OP("+",mod(i+2,L))]))
#         push!(opsList, (Delta*1.0, [OP("Z",i),OP("Z",mod(i+2,L))]))
#         push!(opsList, (-1.0*B_scale, [OP("X",mod(i+1,L))]))
#         push!(opsList, (-g_tau*B_scale, [OP("Z",mod(i+1,L)),OP("Z",mod(i+3,L))]))
#         push!(opsList, (-u_tau*B_scale, [OP("X",mod(i+1,L)),OP("X",mod(i+3,L))]))
#         push!(opsList, (alpha*1.0, [OP("Z",i),OP("X",mod(i+1,L)), OP("Z",mod(i+2,L))]))
#     end
#     return opsList
# end


# if testing2
#     println("Testing abstract basis for XXZ*")
#     L=4
#     Delta = 0.8
#     alpha = 0.0
#     g_tau  = 0.3
#     u_tau = 0.1
#     B_scale = 1.0
#     H_XXZ_1 = make_XXZ_star_operators(L,Delta,alpha,g_tau,u_tau,B_scale)

# end



# function make_XXZ_star_operators_new(L :: Int, Delta :: Float64, alpha :: Float64, g_tau :: Float64, u_tau :: Float64, B_scale :: Float64)

# 	H = HAMILTONIAN("XXZ_star",[])
#     for i in 0:2:L-1
# 		H += term(0.5, 				OP("+", i), 	OP("-", (i+2) % L))
# 		H += term(0.5, 				OP("-", i), 	OP("+", (i+2) % L))
# 		H += term(Delta, 			OP("Z", i), 	OP("Z", (i+2) % L))
# 	 	H += term(-1.0*B_scale, 	OP("X",(i+1) % L))
# 	 	H += term(-g_tau*B_scale, 	OP("Z",(i+1) % L),OP("Z",(i+3) % L))
# 	 	H += term(-u_tau*B_scale, 	OP("X",(i+1) % L),OP("X",(i+3) % L))
# 	 	H += term(alpha, 			OP("Z",i),		OP("X",(i+1) % L),	OP("Z", (i+2)%L))
#     end
#     return H
# end

# if testing2
#     println("Testing abstract basis for XXZ*")
#     L=4
#     Delta = 0.8
#     alpha = 0.0
#     g_tau  = 0.3
#     u_tau = 0.1
#     B_scale = 1.0
#     H_XXZ_2 = make_XXZ_star_operators_new(L,Delta,alpha,g_tau,u_tau,B_scale)
#     [println(TERM(H_XXZ_1[i][1],H_XXZ_1[i][2]),"\t", H_XXZ_2.terms[i]) for i in 9:10]
# end


################ Hamiltonian constructor for a given basis ##########################




function make_Hamiltonian(L :: Int, basis :: Basis, abstract_Ham :: HAMILTONIAN)
    #L - length of the spin chain
    #basis - basis Dict
    #reps - represenatives of conjugacy classes Dict
    #opsList - list of operators we want to add to the Hamiltonian

    dim = length(basis.conj_classes)
    H = spzeros(ComplexF64,dim,dim)

    #iterate over conjugacy classes
    for (x,cc_x) in basis.conj_classes

        #loop over terms in the Hamiltonian
        for term in abstract_Ham.terms
            #find |y> = H |x> ####ASSUMPTION: this is a state, not superposition
            y = apply_operators(VEC(1.0,x),term)

            #check if the element exists --- it could have zero norm
            if ~is_Null_VEC(y) && haskey(basis.get_conj_class,y.state)

                bv_y = basis.get_conj_class[y.state]
                cc_y = basis.conj_classes[bv_y.conj_class]

                #check if upper diagonal
                if cc_y.index >= cc_x.index
                    #compute the matrix element --- see Note
                    h_xy = numchop(sqrt(cc_y.norm/cc_x.norm) * bv_y.phase_factor * y.factor)
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
    #maybe a bug?
    #return Hermitian(H,:U)
end 

if testing2
    println("Testing making the Hamiltonian from the full basis")
    abstract_hamiltonian = make_XXZ_operators_new(4,0.8)
   	#println(abstract_hamiltonian)
    H2 = make_Hamiltonian(4,basisFull,abstract_hamiltonian)
    println(H1 == H2)
end
