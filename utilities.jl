
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



########code for checking if two bases are the same

function check_same_basis_debug(basis1,basis2)

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