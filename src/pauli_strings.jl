################
#
# This file implements operators on PAULI strings, with efficient techniques
# for multiplication, taking commutators, and applying PAULI strings to states
#
################

"""
Internal representation of a pauli string as two bitstrings with coefficient i^delta (-1)^epsilon.
"""
struct VERT
	delta :: Bool
	epsilon :: Bool
	v :: UInt64
	w :: UInt64
end


const zeroVERT = VERT(false,false,UInt64(0),UInt64(0))


# """
# Finds the first non-I place.
# """
function edges(V :: VERT) :: Tuple{Int64,Int64}
	ones = V.v | V.w
	return (ones == 0) ? (0,0) : (trailing_zeros(ones),63-leading_zeros(ones))
end

function shift(V :: VERT) :: VERT
	start = trailing_zeros(V.v | V.w)
	return start == 0 ? V : VERT(V.delta, V.epsilon, V.v >>> start, V.w >>> start)
end

function shift_left(V :: VERT, n :: Int64) :: VERT
    return VERT(V.delta, V.epsilon, V.v >>> n, V.w >>> n)
end

function shift_right(V :: VERT, n :: Int64) :: VERT
    return VERT(V.delta, V.epsilon, V.v << n, V.w << n)
end


function vertString(V :: VERT) :: Tuple{String,String}
	(start, stop) = edges(V)
	numYs = count_ones(V.v & V.w) #count number of (-i)'s out front
	# println(numYs)
	if V.delta
		numYs = numYs +1
	end
	if V.epsilon
		numYs = numYs +2
	end
	numYs = mod(numYs,4)
	# println(numYs)

	vs = digits(V.v,base=2,pad=stop+1)
	ws = digits(V.w,base=2,pad=stop+1)

	names = ["I" "X"; "Z" "Y"]
	arr = reverse([names[vs[i]+1,ws[i]+1] for i in start+1:stop+1])
	opnames = foldr(*,arr)
	leading = start > 0 ? (start > 1 ? "I^$(start)" : "I") : ""

	coeff = ""
	if numYs == 0 
		coeff = ""
	elseif numYs == 1
		coeff = "i*"
	elseif numYs == 2
		coeff = "-1*"
	else
		coeff = "-i*"
	end

	return (coeff,"$(opnames)$(leading)")
end

"""
Writes the string in a human-readable way
"""
function Base.show(io :: IO, V :: VERT)
	if (V.v | V.w) == 0
		print(io,"I")
	else
		# print(io,vertString(V))
		(coeff,verts) = vertString(V)
		print(io,coeff*verts)
	end
end

"""
External representation of PAULI string.
"""
struct PAULI
	x :: Complex{Float64}
	V :: VERT

	function PAULI(name :: String)
		names = reverse(map(String,split(name,"")))
		numYs = 0
		v :: UInt64  = 0
		w :: UInt64  = 0

		for (i,name) in enumerate(names)
			if name == "I"
				(vi, wi) = (0,0)
			else
				if name == "X"
					(vi, wi) = (0,1)
				elseif name == "Y"
					(vi, wi) = (1,1)
					numYs += 1
				elseif name == "Z" 
					(vi, wi) = (1,0)
				else 
					error("invalid name!")
				end
			end
			# println("i: ", i, " vi: ", vi," wi: ", wi)
			v = v | (vi << (i-1)) 
			w = w | (wi << (i-1))
			# println("v: ", v," w:", w)
		end

		numYs = mod(numYs,4)
		delta = (numYs == 1 || numYs == 3)
		epsilon = (numYs == 1 || numYs == 2)


		return new(1, VERT(delta,epsilon,v,w))
	end

	function PAULI(x :: Number, V::VERT)
		return new(x,V)
	end
end


const zeroPAULI = PAULI(1,zeroVERT)


function cannonicalForm(P :: PAULI) :: PAULI
	x = P.x
	if P.V.epsilon
		x = -x
	end
	if P.V.delta
		x = Complex(-x.im,x.re)
	end

	return PAULI(x,VERT(false,false,P.V.v,P.V.w))
end

function shift_left(P :: PAULI, n::Int64)::PAULI PAULI(P.x,shift_left(P.V,n)) end
function shift_right(P :: PAULI, n::Int64)::PAULI PAULI(P.x,shift_right(P.V,n)) end




function Base.:*(z :: Number, P :: PAULI)
	return PAULI(z*P.x,P.V)
end


"""
isNonTrivial(P)

Tests if a PAULIOp is the trivial one,i.e. I_0 I_1...I_63
"""
function isNonTrivial(P :: PAULI) :: Bool
	return P.V.v | P.V.w != 0
end

function string_length(P::PAULI)
	(s,e) =  edges(P.V)
	return e-s+1
end

#check coefficients first, since 0*vert = 0*vert'
function isApprox(P1 :: PAULI, P2 :: PAULI) :: Bool
	P1 = cannonicalForm(P1)
	P2 = cannonicalForm(P2)
	return isapprox(P1.x, P2.x) && P1.V.v == P2.V.v && P1.V.w == P2.V.w
end



function showVerbose(V :: VERT)
	(start,stop) = edges(V)
	print("VERT[$(V.epsilon), $(V.delta), $(digits(V.v,base=2,pad=stop+1)), $(digits(V.w,base=2,pad=stop+1)), $(start), $(stop)]\n")
end

function showVerbose(P :: PAULI)
	(start,stop) = edges(P.V)
	print("PAULI[$(P.x), $(digits(P.V.v,base=2,pad=stop+1)), $(digits(P.V.w,base=2,pad=stop+1)), $(start), $(stop)]\n")
end


"""
Efficient multiplication of PAULI strings
"""
function Base.:*(V1 :: VERT, V2 :: VERT) :: VERT

	d1 = (count_ones(V1.w & V2.v) % 2) == 1
	# d2 = count_ones(V2.w & V1.v) % 2

	#xor = binary addition mod 2
	#and = binary multiplication mod 2
	delta_new = xor(V1.delta, V2.delta)
    d1d2 = (V1.delta && V2.delta)
    
	epsilon_new = xor(V1.epsilon,V2.epsilon,d1d2,d1)
	v_new = xor(V1.v,V2.v)
	w_new = xor(V1.w,V2.w)

	return VERT(delta_new,epsilon_new,v_new,w_new)
end

"""
Efficient multiplication of PAULI strings assuming delta = epsilon = false
"""
function multFast(V1 :: VERT, V2 :: VERT) :: VERT

	d1 = (count_ones(V1.w & V2.v) % 2) == 1

	v_new = xor(V1.v,V2.v)
	w_new = xor(V1.w,V2.w)

	return VERT(false,d1,v_new,w_new)
end


"""

Note: if [alpha, beta] commutes, then returns I, since 0 doesn't make sense as a PAULIInt
When PAULIInt is converted to PAULI, this becomes the zero PAULI
"""
function commutator(V1 :: VERT, V2 :: VERT) :: VERT


	# d1 = (count_ones(V1.w & V2.v) % 2) == 1
	# d2 = (count_ones(V2.w & V1.v) % 2) == 1
	# println(d1)
	# println(d2)

	D1 = V1.w & V2.v
	d = count_ones(xor(D1, V2.w & V1.v)) % 2 == 1
	# println(d)

	# println("V1: ", V1, " V2: ", V2)
	if d
		d1 = (count_ones(D1) % 2) == 1

		#xor = binary addition mod 2
		#and = binary multiplication mod 2
		delta_new = xor(V1.delta, V2.delta)
	    d1d2 = (V1.delta && V2.delta)
	    
		epsilon_new = xor(V1.epsilon,V2.epsilon,d1d2,d1)
		v_new = xor(V1.v,V2.v)
		w_new = xor(V1.w,V2.w)

		return VERT(delta_new,epsilon_new,v_new,w_new)
	else
		return zeroVERT
	end
end

const PauliOp = Dict{VERT,Complex{Float64}}

function cannonicalForm(x :: Complex{Float64}, V :: VERT) :: Tuple{Complex{Float64},VERT}
	if V.epsilon
		x = -x
	end
	if V.delta
		x = Complex(-x.im,x.re)
	end

	return (x,VERT(false,false,V.v,V.w))
end


function cannonicalForm(O :: PauliOp) :: PauliOp
	O_new = Dict{VERT,Complex{Float64}}()
	for (V,x) in O
		(x,V) = cannonicalForm(x,V)
		O[V] = haskey(O,V) ? O[V] + x : x
	end

	return O_new
end


function Base.:*(O1 :: PauliOp, O2 :: PauliOp) :: PauliOp
	
	O = Dict{VERT,Complex{Float64}}()

	for (V1,x1) in O1, (V2,x2) in O2
		(x,V) = cannonicalForm(x1*x2,V1*V2)
		O[V] = haskey(O,V) ? O[V] + x : x
	end

	return O

end

function Base.:+(O1 :: PauliOp, O2 :: PauliOp) :: PauliOp
	# O1 = cannonicalForm(O1)
	# O2 = cannonicalForm(O2)

	# return merge(+,O1,O2) #for some reason this
	#doesn't copy the dict correctly. doing it manually instead...

	O = PauliOp()

	for (V1,x1) in O1
		(x,V) = cannonicalForm(x1,V1)
		O[V] = haskey(O,V) ? O[V] + x : x
	end

	for (V2,x2) in O2
		(x,V) = cannonicalForm(x2,V2)
		O[V] = haskey(O,V) ? O[V] + x : x
	end
	return O

end

function Base.:*(z :: Number, O :: PauliOp) :: PauliOp
	return Dict( V =>  z*x for (V,x) in O)
end

function string_length(O :: PauliOp)
	return maximum( string_length(V) for (V,x) in O)
end

function adj(x :: Complex{Float64}, V :: VERT) :: Pair{VERT,Complex{Float64}}
	x = conj(x)
	flip_ones =  count_ones(V.v & V.w)  & 1 == 1
	epsilon = xor(V.delta,V.epsilon,flip_ones)
	return VERT(V.delta,epsilon,V.v,V.w) => x
end

function adj(O :: PauliOp) :: PauliOp
	return Dict( adj(x,V) for (V,x) in O)
end

function trim(O :: PauliOp; EPS = 10e-15) :: PauliOp
	return filter(Vx->abs(Vx[2]) > EPS, O)
end




# function commutator(V1 :: VERT, V2 :: VERT, offset :: Int64) :: PAULIInt
# 	# println("offset: ", offset)
# 	if offset < 0
# 		V2 = VERT(V2.v >>> offset, V2.w >>> offset)
# 	elseif offset > 0
# 		V1 = VERT(V1.v << offset, V1.w << offset)
# 	end

# 	# println("V1: ",vertString(V1))
# 	# println("V2: ",vertString(V2))


# 	d1 = count_ones(V1.w & V2.v) % 2
# 	d2 = count_ones(V2.w & V1.v) % 2

# 	# println("V1: ", V1, " V2: ", V2)
# 	if d1 != d2
# 		epsilon = (d1 == 1)
# 		new_v = xor(V1.v,V2.v)
# 		new_w = xor(V1.w,V2.w)

# 		#shift over into cannonical form.
# 		#only makes sense with translation invariance
# 		# println(VERT(new_v,new_w))
# 		V = shift(VERT(new_v,new_w))
# 		p = PAULIInt(false,epsilon,V)
# 		# println(p)
# 		return p
# 	else
# 		return zeroPAULI
# 	end

# end


# function convolve(P1 :: PAULI, P2 :: PAULI) :: Op

# 	V1 = P1.V
# 	V2 = P2.V

# 	e1 = edges(V1)
# 	e2 = edges(V2)

# 	min_offset = e2[1] - e1[2]
# 	max_offset = e2[2] - e1[1]

# 	pre = P1.x * P2.x

# 	# return map(P-> PAULI(pre,P),
# 	# 		 filter(isNonTrivial,
# 	# 		  [commutator(V1,V2,offset)
# 	# 		   for offset in min_offset:max_offset]
# 	# 		   	)
# 	# 		 )

# 	return [PAULI(pre,P) for P in 
# 			[commutator(V1,V2,offset) for offset in min_offset:max_offset]
# 			if isNonTrivial(P)]

# 	# return [commutator(V1,V2,offset)
# 			   # for offset in min_offset:max_offset]
# end



# function OP(d :: Dict{String,T}) :: Op where {T<:Number}
# 	return [coefficient*PAULI(name) for (name,coefficient) in d]
# end



# function L(P :: PAULI, H :: Op) :: Op
# 	# println("P: ", P)
# 	# println("H: ", H)

# 	newTerms = Op([])
# 	# terms = Dict{VERT,PAULI}()

# 	for h in H
# 		# println("current term in H: ", h)
# 		# c = convolve(h,P)
# 		# println(c)
# 		newTerms = [newTerms; convolve(h,P)]
# 		# c = convolve(h,P)
# 		# for p in c
# 		# 	v = p.V
# 		# 	if !haskey(terms,v)
# 		# 		terms[v] = p 
# 		# 	else 
# 		# 		p0 = terms[v]
# 		# 		pre_new = p0.x + p.x
# 		# 		terms[v] = PAULI(pre_new,v)
# 		# 	end
# 		# end
# 	end
# 	# println("[H,P]:", newTerms)
# 	# return collect(values(terms))
# 	return newTerms
# end

# function T(P :: PAULI, H :: Op) :: Set{VERT}
# 	# pr intln("P: ", P)
# 	# println("H: ", H)

# 	# newTerms = Op([])
# 	terms = Dict{VERT,Complex{Float64}}()
# 	# s = Set{VERT}()
# 	for h in H
# 		# println("current term in H: ", h)
# 		c = convolve(h,P)
# 		# println(c)
# 		# newTerms = [newTerms; convolve(h,P)]
# 		# c = Set(p.V for p in convolve(h,P))
# 		# s = union(c,s)
# 		for p in c
# 			v = p.V
# 			if !haskey(terms,v)
# 				terms[v] = p.x 
# 			else 
# 				p0 = terms[v]
# 				pre_new = p0 + p.x
# 				terms[v] = pre_new
# 			end
# 		end
# 	end
# 	# println("[H,P]:", newTerms)
# 	return Set(vert for (vert, coeff) in terms if !(abs(coeff)<10e-12))
# 	# return s
# end

# function VEC(O :: Op) :: Vec
# 	return Dict(P.V => P.x for P in O)
# end

# function vec_norm(v :: Vec) :: Float64
# 	if length(v) == 0
# 		return 0.0
# 	else 
# 		return sqrt(sum(abs2(val) for val in values(v)))
# 	end
# end

# function vec_scalar(x :: T, v :: Vec) :: Vec where {T<:Number}
# 	return Dict(vert => x*coeff for (vert,coeff) in v)
# end

# #v1^dag .v2
# function vec_ip(v1 :: Vec, v2 :: Vec) :: Complex{Float64}
# 	n = 0.0
# 	for (vert, c1) in v1
# 		if haskey(v2,vert)
# 			c2 = v2[vert]
# 			n += conj(c1)*c2
# 		end
# 	end
# 	return n
# end

# function vec_subtract(v :: Vec, w :: Vec) :: Vec
# 	z = copy(v)
# 	for (key_w,val_w) in w
# 		if !haskey(z,key_w)
# 			z[key_w] = -val_w
# 		else
# 			z[key_w] -= val_w
# 		end	
# 	end
# 	return z
# end

# function vec_subtract_scalar(v :: Vec, lambda :: T, w :: Vec) :: Vec where {T <: Number}
# 	z = copy(v)
# 	for (key_w,val_w) in w
# 		if !haskey(z,key_w)
# 			z[key_w] = -lambda*val_w
# 		else
# 			z[key_w] -= lambda*val_w
# 		end	
# 	end
# 	return z
# end

# function multMatrixFree(H :: Op, v :: Vec) :: Vec
# 	w = Vec()
# 	for (vert, coeff) in v
# 		if abs(coeff) >= 10^(-12)
# 			P = PAULI(coeff,vert)
# 			hP = L(P,H)

# 			for P_new in hP
# 				new_vert = P_new.V
# 				new_coeff = P_new.x

# 				if !haskey(w,new_vert)
# 					w[new_vert] = new_coeff
# 				else
# 					w[new_vert] += new_coeff
# 				end
# 				# println(new_vert, ": ", new_coeff)
# 			end
# 		end	
# 	end
# 	return w
# end

# function multMatrixFreeSemiclassical(H :: Op, v :: Vec) :: Vec
# 	w = Vec()
# 	for (vert, coeff) in v
# 		# println(coeff)
# 		if abs(coeff) >= 10^(-14)
# 			P = PAULI(coeff,vert)
# 			hP = T(P,H)
# 			# println("T($(P), $(H)) = ", hP)

# 			for v_new in hP
# 				# println(v_new,": ", coeff)
# 				if !haskey(w,v_new)
# 					w[v_new] = coeff
# 				else
# 					w[v_new] += coeff
# 				end
# 			end
# 		end
# 	end
# 	return w
# end

# function LanczosMatrixFree(H :: Op, O_0 :: Op, m :: Int) :: Array{Float64}

# 	betas = Array{Float64}(undef,m)

# 	v = VEC(O_0)
# 	n0 = vec_norm(v)
# 	v = vec_scalar(1/n0,v)

# 	v_old = copy(v)
# 	w = multMatrixFree(H,v)
# 	println("Vector 1: size $(length(w))")
# 	# println(w)
# 	beta = vec_norm(w)
# 	betas[1] = beta

# 	for j in 2:m
# 		v_old = vec_scalar(beta,v)
# 		v = vec_scalar(1/beta,w) #v_j = w_{j-1}/beta_{j-1}
# 		w = multMatrixFree(H,v)
# 		w = vec_subtract(w,v_old)
# 		println("Vector $(j): size $(length(w))")
# 		# println(w)
# 		beta = vec_norm(w)
# 		betas[j] = beta
# 	end

# 	return betas
# end


# function LanczosMatrixFreeVectors(
# 		H :: Op, 
# 		O_0 :: Op,
# 		m :: Int
# 		) :: Tuple{Array{Float64},Array{Vec}}

# 	betas = Array{Float64}(undef,m)
# 	vs = Array{Vec}(undef,m)

# 	v = VEC(O_0)
# 	n0 = vec_norm(v)
# 	v = vec_scalar(1/n0,v)
# 	vs[1] = copy(v)
# 	println("Starting vec:", v)



# 	v_old = copy(v)
# 	w = multMatrixFree(H,v)
# 	println("Vector 1: size $(length(w))")
# 	# println(w)
# 	beta = vec_norm(w)
# 	betas[1] = beta

# 	for j in 2:m
# 		v_old = vec_scalar(beta,v)
# 		v = vec_scalar(1/beta,w) #v_j = w_{j-1}/beta_{j-1}
# 		vs[j] = copy(v)
# 		w = multMatrixFree(H,v)


# 		# println("IP (j+1,j):", vec_ip(w,v))
# 		w = vec_subtract(w,v_old)
# 		println("Vector $(j): size $(length(w))")
# 		# println(w)
# 		beta = vec_norm(w)
# 		betas[j] = beta
# 	end

# 	return (betas,vs)
# end


# function LanczosMatrixFreeVectorsSemiclassical(
# 		H :: Op, 
# 		O_0 :: Op,
# 		m :: Int
# 		) :: Tuple{Array{Complex{Float64}},Array{Float64},Array{Vec}}

# 	alphas = Array{Complex{Float64}}(undef,m+1)
# 	betas = Array{Float64}(undef,m)
# 	vs = Array{Vec}(undef,m)

# 	v = VEC(O_0)
# 	#make it all real, otherwise the Y's coefficients are wrong
# 	v = Vec(vert => Complex{Float64}(abs(coeff)) for (vert, coeff) in v)
# 	n0 = vec_norm(v)
# 	v = vec_scalar(1/n0,v)
# 	vs[1] = copy(v)
# 	# println("starting vec:", v)
# 	# println("values: ", values(v))

# 	v_old = copy(v)
# 	w = multMatrixFreeSemiclassical(H,v)
# 	# println("Vector 1: size $(length(w))")
# 	# println(w)

# 	alpha = vec_ip(w,v)
# 	# println("alpha: (2,1):", alpha)
# 	alphas[1] = alpha
# 	# v = vec_scalar(alpha,v)
# 	w = vec_subtract_scalar(w,alpha,v)

# 	beta = vec_norm(w)
# 	betas[1] = beta

# 	for j in 2:m
# 		v_old = v#vec_scalar(beta,v) #v_{j-1}

# 		v = vec_scalar(1/beta,w) #v_j = w_{j-1}/beta_{j-1}
# 		vs[j] = copy(v)
# 		w = multMatrixFreeSemiclassical(H,v)

# 		alpha = vec_ip(w,v)
# 		# println("alpha: ($(j+1), $(j)) ", alpha)
# 		alphas[j+1] = alpha
# 		# v = vec_scalar(alpha,v)
# 		# w = vec_subtract(w,v)
# 		# v = vec_scalar(1/alpha,v)
# 		# w = vec_subtract(w,v_old)
# 		w = vec_subtract_scalar(w, alpha, v)
# 		w = vec_subtract_scalar(w, beta, v_old)

# 		println("Vector $(j): size $(length(w))")
# 		# println(w)
# 		beta = vec_norm(w)
# 		betas[j] = beta
# 	end
# 	return (alphas,betas,vs)
# end

# """
# Calculates the length of the support needed for the vertex.
# """
# function support(V :: VERT) :: Int64
# 	(start, stop) = edges(V)
# 	return stop-start
# end

# """
# Calculates the number of non-I operators, i.e. the size of the operator.
# """
# function vert_size(V :: VERT) :: Int64
# 	ones = (V.v | V.w)
# 	s = count_ones(ones)
# 	return s
# end


# """
# Counts the lengths of a vector. The length is ``l=end - begin``.
# """
# function support_dist(v :: Vec) :: Array{Int64}
# 	ls = zeros(Int64,64)

# 	for (vert,coeff) in v
# 		if !isapprox(coeff,0.0+0.0im)
# 			l = support(vert)
# 			ls[l+1] += 1
# 		end
# 	end

# 	return ls
# end

# """
# Counts the lengths of a vector. The length is ``l=end - begin``.
# """
# function size_dist(v :: Vec) :: Array{Int64}
# 	ss = zeros(Int64,64)

# 	for (vert,coeff) in v
# 		if !isapprox(coeff,0.0+0.0im)
# 			s = vert_size(vert)
# 			ss[s] += 1
# 		end
# 	end

# 	return ss
# end

# """
# Counts the lengths of a vector. The length is ``l=end - begin``.
# """
# function support_dist_weighted(v :: Vec) :: Array{Float64}
# 	ls = zeros(Float64,64)

# 	for (vert,coeff) in v
# 		if !isapprox(coeff,0.0+0.0im)
# 			l = support(vert)
# 			ls[l+1] += abs2(coeff)
# 		end
# 	end

# 	return ls
# end

# """
# Counts the lengths of a vector. The length is ``l=end - begin``.
# """
# function size_dist_weighted(v :: Vec) :: Array{Float64}
# 	ss = zeros(Float64,64)

# 	for (vert,coeff) in v
# 		if !isapprox(coeff,0.0+0.0im)
# 			s = vert_size(vert)
# 			ss[s] += abs2(coeff)
# 		end
# 	end

# 	return ss
# end

