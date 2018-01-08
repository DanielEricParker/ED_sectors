
#####################################################################
#Miscellaneous utilities for pretty printing, and little tasks

#####################################################################



#################### Display/pretty printing ##########################

"""
Pretty printing for spin-1/2 states.
#Arguments
* 'psi :: Array{UInt64}': states to display
"""
function display_states(
	s::Array{UInt64}
	)
    println(map(bin,s))
end


################### General little useful functions #############


# function numchop(x :: Complex{Float64}, epsilon = 10e-13)
#     #x - complex number
#     #returns - x with small real and imaginary parts removed
#     a = real(x)
#     b = imag(x)

#     w = abs(a) < epsilon ? 0 : a
#     z = abs(b) < epsilon ? 0 : b 

#     return complex(w,z)
# end

#borrowed from https://github.com/jlapeyre/ZChop.jl/blob/master/src/ZChop.jl

const zeps = 1e-14

zchop(x::T, eps=zeps) where T <: Real = abs(x) > convert(T,eps) ? x : zero(T)
zchop(x::T, eps=zeps) where T<:Integer = abs(x) > eps ? x : zero(T)
zchop(x::T, eps=zeps) where T<:Complex = complex(zchop(real(x),eps),zchop(imag(x),eps))
zchop(a::T, eps=zeps) where T<:AbstractArray = (b = similar(a); for i in 1:length(a) b[i] = zchop(a[i],eps) end ; b)
zchop!(a::T, eps=zeps)  where T<:AbstractArray = (for i in 1:length(a) a[i] = zchop(a[i],eps) end ; a)
zchop(x::T,eps=zeps) where T<:Union{AbstractString,Char}=  x
zchop(x::T,eps=zeps) where {T<:Irrational} = zchop(float(x),eps)
zchop(x::Expr,eps=zeps) = Expr(x.head,zchop(x.args)...)
zchop(x,eps) =  applicable(start,x) ? map((x)->zchop(x,eps),x) : x
zchop(x) = applicable(start,x) ? map(zchop,x) : x

############### Exporting and File IO ######################

"""
Exports arrays to csv files, with parameters in the filename. Particularly useful when varying 1 or more parameters.
#Arguments
* 'data :: Array{Float64, N} where N': data to export, an Array
* 'name :: String': name/prefix for the data, e.g. "XXZ_eigenvalues"
* 'parameters :: Array{Tuple{String,Number}}': Array of Tuples (parameter, value), e.g. ("Delta", 0.6)
"""
function export_data(
	data :: Array{Float64, N} where N,
	name :: String,
	parameters :: Array{Tuple{String,Number}})
    
    fileName = name
    for (param,value) in parameters
        if typeof(value) == Int
            fileName *= string("_", param, "_", value)
        elseif typeof(value) == Float64 
             fileName *= @sprintf "_%s_%.4f" param value
        end
    end
    fileName *= ".csv"

    #export as CSV
    writedlm(fileName,data,", ")
end



"""
Checks if two bases are the same. Useful for debugging.
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