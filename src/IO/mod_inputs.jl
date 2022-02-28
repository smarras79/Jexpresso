using Crayons.Box
using Revise

export mod_inputs_user_inputs
export mod_inputs_print
export mod_inputs_print_welcome

include("./user_inputs.jl")

function mod_inputs_user_inputs()#(inputs::AbstractDict, nvars::TInt)
    
    inputs = user_inputs() # user_inputs is a Dict
    
    #
    # Correct quantities based on a hierarchy of input variables
    #
    if (inputs[:nsd] == 1)
        inputs[:npy] = 1
        inputs[:npz] = 1
    elseif(inputs[:nsd] == 2)
        inputs[:npz] = 1
    end
    
    nvars::Int8 = 1
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
    else
        println( "\n !!!! ERROR in user_inputs.jl: equation_set ", inputs[:equation_set], " is not coded!!!")
        println( " !!!! Chose among: \n !!!!  [1] burgers1 \n")
        error()
    end
    
    return inputs, nvars
end


function mod_inputs_print(inputs::Dict{}; nvars::Int8)
    
    print(GREEN_FG(" # User inputs:\n"))
    println( " # Equation set:     ", inputs[:equation_set])
    println( " # Problem:          ", inputs[:problem])
    println( " # N. variables:     ", nvars)
    println( " # N. space dims:    ", inputs[:nsd])
    println( " # Polynomial order: ", inputs[:nop])
    println( " # N. x-points:      ", inputs[:npx])
    println( " # [xmin, xmax]:     ", inputs[:xmin], " ", inputs[:xmax])
    if (inputs[:nsd] > 1)
        println( " # N. y-points:      ", inputs[:npy])
        println( " # [ymin, ymax]:     ", inputs[:ymin], " ", inputs[:ymax])
    end
    if (inputs[:nsd] == 3)
        println( " # N. z-points:      ", inputs[:npz])
        println( " # [zmin, zmax]:     ", inputs[:zmin], " ", inputs[:zmax])
    end
    
end


function mod_inputs_print_welcome()

    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))
    print(BLUE_FG(" # Welcome to ", RED_FG("jexpresso\n")))
    print(BLUE_FG(" # A Julia code to solve turbulence problems in the atmosphere\n"))
    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))

end


inputs        = Dict{}()
inputs, nvars = mod_inputs_user_inputs()

