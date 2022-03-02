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
    
    print(GREEN_FG(" # User inputs: ........................ \n"))

    #=for key in sort(collect(keys(inputs)))
        println("$key => $(inputs[key])")
    end=#
    
    
    if (haskey(inputs, :equation_set)) println( " # Equation set:     ", inputs[:equation_set]) end
    
    if (haskey(inputs, :problem))    println( " # Problem:          ", inputs[:problem]) end
    if (haskey(inputs, :lread_gmsh)) println( " # Read GMSH:        ", inputs[:lread_gmsh]) end
    if (haskey(inputs, :nsd))        println( " # N. space dims:    ", inputs[:nsd]) end
    if (haskey(inputs, :nvars))      println( " # N. variables:     ", nvars) end
    if (haskey(inputs, :nop))        println( " # Polynomial order: ", inputs[:nop]) end
    if (haskey(inputs, :npx))        println( " # N. x-points:      ", inputs[:npx]) end
    if (haskey(inputs, :xmin) && haskey(inputs, :xmax)) println( " # [xmin, xmax]:     ", inputs[:xmin], " ", inputs[:xmax]) end
    if ((haskey(inputs, :nsd)) && inputs[:nsd] > 1)
        if (haskey(inputs, :npy)) println( " # N. y-points:      ", inputs[:npy]) end
        if (haskey(inputs, :ymin) && haskey(inputs, :ymax)) println( " # [ymin, ymax]:     ", inputs[:ymin], " ", inputs[:ymax]) end
    end
    if ((haskey(inputs, :nsd)) && inputs[:nsd] == 3)
        if (haskey(inputs, :npz)) println( " # N. z-points:      ", inputs[:npz]) end
        if (haskey(inputs, :zmin) && haskey(inputs, :zmax)) println( " # [zmin, zmax]:     ", inputs[:zmin], " ", inputs[:zmax]) end
    end
     print(GREEN_FG(" # User inputs: ........................ DONE\n"))
end


function mod_inputs_print_welcome()

    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))
    print(BLUE_FG(" # Welcome to ", RED_FG("jexpresso\n")))
    print(BLUE_FG(" # A Julia code to solve turbulence problems in the atmosphere\n"))
    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))

end


inputs        = Dict{}()
inputs, nvars = mod_inputs_user_inputs()

