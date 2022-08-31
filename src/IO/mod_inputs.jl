using Crayons.Box
using PrettyTables
using Revise

export mod_inputs_user_inputs
export mod_inputs_print_welcome

include("./user_inputs.jl")

function mod_inputs_user_inputs()

    error_flag::Int8 = 0
    
    inputs = user_inputs() # user_inputs is a Dict

    print(GREEN_FG(" # User inputs from ...IO/user_inputs.jl .............. \n"))
    pretty_table(inputs; sortkeys=true, border_crayon = crayon"yellow")    
    print(GREEN_FG(" # User inputs: ................................... DONE\n"))    
    #
    # Check that necessary inputs exist in the Dict inside .../IO/user_inputs.jl
    #
    mod_inputs_check(inputs, :equation_set, "e")
    mod_inputs_check(inputs, :problem, "e")
    mod_inputs_check(inputs, :nop, Int8(4), "w") #Polynomial order

    if(!haskey(inputs, :lexact_integration))
        inputs[:lexact_integration] = true #Default integration rule is exact
    end
    
    #Grid entries:
    if(!haskey(inputs, :lread_gmsh) || inputs[:lread_gmsh] == false)
        
        mod_inputs_check(inputs, :nsd,  "e")
        mod_inputs_check(inputs, :nelx,  "e")
        mod_inputs_check(inputs, :xmin, "e")
        mod_inputs_check(inputs, :xmax, "e")
        mod_inputs_check(inputs, :nely,  "e")
        mod_inputs_check(inputs, :ymin, "e")
        mod_inputs_check(inputs, :ymax, "e")
        mod_inputs_check(inputs, :nelz,  "e")
        mod_inputs_check(inputs, :zmin, "e")
        mod_inputs_check(inputs, :zmax, "e")
        
    else
        mod_inputs_check(inputs, :gmsh_filename, "e")
        
        mod_inputs_check(inputs, :nsd,  Int8(3), "-")
        mod_inputs_check(inputs, :nelx,  Int8(2), "-")
        mod_inputs_check(inputs, :xmin, Float64(-1.0), "-")
        mod_inputs_check(inputs, :xmax, Float64(+1.0), "-")
        mod_inputs_check(inputs, :nely,  Int8(2), "-")
        mod_inputs_check(inputs, :ymin, Float64(-1.0), "-")
        mod_inputs_check(inputs, :ymax, Float64(+1.0), "-")
        mod_inputs_check(inputs, :nelz,  Int8(2), "-")
        mod_inputs_check(inputs, :zmin, Float64(-1.0), "-")
        mod_inputs_check(inputs, :zmax, Float64(+1.0), "-")

        s= """ 
           jexpresso: Some undefined (but unnecessary) user inputs 
           MAY have been given some default values.
           User needs not to worry about them.
           """
        @warn s
        
    end #lread_gmsh =#
    
    
    #
    # Correct quantities based on a hierarchy of input variables
    #
    # Define default npx,y,z for native grid given
    # values for the user's nelx,y,z
    if(haskey(inputs, :nelx))
        inputs[:npx] = inputs[:nelx] + 1
    else
        inputs[:npx] = 2
    end
    if(haskey(inputs, :nely))
        inputs[:npy] = inputs[:nely] + 1
    else
        inputs[:npy] = 2
    end
    if(haskey(inputs, :nelz))
        inputs[:npz] = inputs[:nelz] + 1
    else
        inputs[:npz] = 2
    end
    
    if (inputs[:nsd] == 1)
        inputs[:npy] = 1
        inputs[:npz] = 1
    elseif(inputs[:nsd] == 2)
        inputs[:npz] = 1
    end

    #
    # Define nvars based on the problem being solved
    #
    nvars::Int8 = 1
    if (lowercase(inputs[:equation_set]) == "burgers")
        if(inputs[:nsd] == 1)
            nvars = 1
        elseif (inputs[:nsd] == 2)
            nvars = 2
        end
        println( " # nvars     ", nvars)
        
    elseif (lowercase(inputs[:equation_set]) == "ns")
        if (inputs[:nsd] == 1)
            nvars = 3
        elseif(inputs[:nsd] == 2)
            nvars = 4
        elseif(inputs[:nsd] == 3)
            nvars == 5
        end
        println( " # nvars     ", nvars)
    else
        s = """
            jexpresso  user_inputs.jl: equation_set ", inputs[:equation_set], " is not coded!
            Chose among:
                    [1] BURGERS
                    [2] NS
        """
        @error s
    end
    
    
    return inputs, nvars
end

function mod_inputs_check(inputs::Dict, key, error_or_warning::String)
    
    if (!haskey(inputs, key))
        s = """
            jexpresso: $key is missing in .../IO/user_inputs.jl
            """
        if (error_or_warning=="e")
            error(s)
        elseif (error_or_warning=="w")
            @warn s
        end
        error_flag = 1
    end
       
end


function mod_inputs_check(inputs::Dict, key, value, error_or_warning::String)

    if (!haskey(inputs, key))
        s = """
            jexpresso: $key is missing in .../IO/user_inputs.jl
            The default value $key=$value will be used.
            """
        if (error_or_warning=="e")
            error(s)
        elseif (error_or_warning=="w")
            @warn s     
        end
        
        #assign a dummy default value
        inputs[key] = value
    end

end

function mod_inputs_print_welcome()

    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))
    print(BLUE_FG(" # Welcome to ", RED_FG("jexpresso\n")))
    print(BLUE_FG(" # A Julia code to solve turbulence problems in the atmosphere\n"))
    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))

end
