
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
    println(map(st -> digits(st,base=2),s))
end


################### General little useful functions #############


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

    export_data(data,name,parameters)

Exports arrays to csv files, with parameters in the filename. Particularly useful when varying 1 or more parameters.
#Arguments
* 'data :: Array{Float64, N} where N': data to export, an Array
* 'name :: String': name/prefix for the data, e.g. "XXZ_eigenvalues"
* 'parameters :: Array{Tuple{String,Number}}': Array of Tuples (parameter, value), e.g. ("Delta", 0.6)
"""
function export_data(
	# data :: Array{Float64, N} where N,
	# name :: String,
	# parameters :: Array{Tuple{String,Number}})
    data, 
    name,
    parameters;
    verbose = false
    )
    fileName = name
    for (param,value) in parameters
        if typeof(value) == Int
            fileName *= string("_", param, "_", value)
        elseif typeof(value) == Float64 
             fileName *= @sprintf "_%s_%.4f" param value
        elseif typeof(value) == String
            fileName *= string("_", param, "_", value)
        end
    end
    fileName *= ".csv"

    if verbose
        println("Exporting data to '$(fileName)'.")
    end
    #export as CSV
    writedlm(fileName,data,", ")

    return fileName
end

"""
    
    make_filename(name,parameters)

Given the prefix 'name' and the parameter/value pairs, generates a filename for the data.

#Arguments
* 'name :: String': name/prefix for the data, e.g. "XXZ_eigenvalues"
* 'parameters :: Array{Tuple{String,Number}}': Array of Tuples (parameter, value), e.g. ("Delta", 0.6)
"""
function make_filename(
    name,
    parameters
    )
    fileName = name
    for (param,value) in parameters
        if typeof(value) == Int
            fileName *= string("_", param, "_", value)
        elseif typeof(value) == Float64 
             fileName *= @sprintf "_%s_%.4f" param value
        elseif typeof(value) == String
            fileName *= string("_", param, "_", value)
        end
    end
    fileName *= ".csv"

    return fileName
end

