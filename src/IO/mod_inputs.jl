include("./user_inputs.jl")

export mod_inputs_user_inputs

function mod_inputs_user_inputs(TInt, TFloat)
    
    inputs = user_inputs()

    #
    # Correct quantities based on a hierarchy of input variables
    #
    if (inputs[:nsd] == 1)
        inputs[:npy] = 1
        inputs[:npz] = 1
    elseif(inputs[:nsd] == 2)
        inputs[:npy] = 1
    end
    
    nvars::TInt = 0
    if (lowercase(inputs[:equation_set]) == "burgers")
        if(inputs[:nsd] == 1)
            nvars = 1
        elseif (inputs[:nsd] == 2)
            nvars = 2
        end
    elseif (lowercase(inputs[:equation_set]) == "ns")
        if (inputs[:nsd] == 1)
            nvars = 3
        elseif(inputs[:nsd] == 2)
            nvars = 4
        elseif(inputs[:nsd] == 3)
            nvars == 5
        end
    end
            
    return inputs, nvars
end
        
