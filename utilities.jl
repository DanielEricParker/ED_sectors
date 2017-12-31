#####################################################################
#Miscellaneous utilities for pretty printing, and little tasks

#####################################################################



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


############### Exporting and File IO ######################

function export_data(data,name,parameters)
    #data - data to export, an Array
    #name - name/prefix for the data, e.g. "XXZ_eigenvalues"
    #parameters - Array of Tuples (parameter, value), e.g. ("Delta", 0.6)
    
    fileName = name
    for (param,value) in parameters
        fileName *= string("_", param, "_", value)
    end

    #export as CSV
    writedlm(filename,data,", ")
end
