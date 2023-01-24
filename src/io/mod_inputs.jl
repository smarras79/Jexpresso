using Crayons.Box
using PrettyTables
using Revise

export mod_inputs_user_inputs
export mod_inputs_print_welcome

include("../../user_inputs.jl")

function mod_inputs_user_inputs()
    
    error_flag::Int8 = 0
    
    inputs = user_inputs() # user_inputs is a Dict
    
    print(GREEN_FG(" # User inputs from ...IO/user_inputs.jl .............. \n"))
    pretty_table(inputs; sortkeys=true, border_crayon = crayon"yellow")    
    print(GREEN_FG(" # User inputs: ................................... DONE\n"))    

    #
    # Check that necessary inputs exist in the Dict inside .../IO/user_inputs.jl
    #
    mod_inputs_check(inputs, :problem, "e")
    mod_inputs_check(inputs, :nop, Int8(4), "w")  #Polynomial order
    
    #Time:
    if(!haskey(inputs, :diagnostics_interval))
        inputs[:diagnostics_interval] = Int8(1)
    end
    mod_inputs_check(inputs, :tend, "e") #Final time
    mod_inputs_check(inputs, :Δt, Float64(1.0), "w") #Δt --> this will be computed from CFL later on
    if(!haskey(inputs, :tinit))
        inputs[:tinit] = 0.0  #Initial time is 0.0 by default
    end
    
    if(!haskey(inputs, :lexact_integration))
        inputs[:lexact_integration] = false #Default integration rule is INEXACT
    end

    mod_inputs_check(inputs, :interpolation_nodes, String("lgl"), "w")
    if(haskey(inputs, :interpolation_nodes) && inputs[:interpolation_nodes] == "llg" || inputs[:interpolation_nodes] == "gll")
        inputs[:interpolation_nodes] = "lgl"
    end
    
    #Grid entries:
    if(!haskey(inputs, :lread_gmsh) || inputs[:lread_gmsh] == false)
        
        mod_inputs_check(inputs, :nsd,  Int8(1), "-")
        mod_inputs_check(inputs, :nelx, "e")
        mod_inputs_check(inputs, :xmin, "e")
        mod_inputs_check(inputs, :xmax, "e")
        mod_inputs_check(inputs, :nely,  Int8(2), "-")
        mod_inputs_check(inputs, :ymin, Float64(-1.0), "-")
        mod_inputs_check(inputs, :ymax, Float64(+1.0), "-")
        mod_inputs_check(inputs, :nelz,  Int8(2), "-")
        mod_inputs_check(inputs, :zmin, Float64(-1.0), "-")
        mod_inputs_check(inputs, :zmax, Float64(+1.0), "-")
        
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
    # Some physical constants and parameters:
    #    
    if(!haskey(inputs, :νx))
        inputs[:νx] = Float16(0.0) #default kinematic viscosity
    end
    if(!haskey(inputs, :νy))
        inputs[:νy] = Float16(0.0) #default kinematic viscosity
    end
    if(!haskey(inputs, :νz))
        inputs[:νz] = Float16(0.0) #default kinematic viscosity
    end
    
    #
    # Correct quantities based on a hierarchy of input variables
    #
    # Define default npx,y,z for native grid given
    # values for the user's nelx,y,z
    if(haskey(inputs, :nelx))
        inputs[:npx] = inputs[:nelx] + 1
    else
        inputs[:npx] = Int8(2)
    end
    if(haskey(inputs, :nely))
        inputs[:npy] = inputs[:nely] + 1
    else
        inputs[:npy] = Int8(2)
    end
    if(haskey(inputs, :nelz))
        inputs[:npz] = inputs[:nelz] + 1
    else
        inputs[:npz] = Int8(2)
    end
    
    if (inputs[:nsd] == 1)
        inputs[:npy] = Int8(1)
        inputs[:npz] = Int8(1)
    elseif(inputs[:nsd] == 2)
        inputs[:npz] = Int8(1)
    end

        
    """
    To add a new set of governing equations, add a new :problem key 
    and add the following lines 

     elseif (lowercase(inputs[:problem]) == "new problem name")
        inputs[:problem] = NewProblemName()
        
        nvars = INTEGER VALUE OF THE NUMBER OF UNKNOWNS for this problem.
        println( " # nvars     ", nvars)

    """
    
    #------------------------------------------------------------------------
    # Define nvars based on the problem being solved
    #------------------------------------------------------------------------
    nvars::Int8 = 1
    if (lowercase(inputs[:problem]) == "burgers")
        inputs[:problem] = burgers()
        
        if(inputs[:nsd] == 1)
            nvars = 1
        elseif (inputs[:nsd] == 2)
            nvars = 2
        end
        inputs[:nvars] = nvars
        println( " # nvars     ", nvars)
        
    elseif (lowercase(inputs[:problem]) == "sw")
        inputs[:problem] = sw()
        
        if (inputs[:nsd] == 1)
            nvars = 2
        elseif(inputs[:nsd] == 2)
            nvars = 3
        elseif(inputs[:nsd] == 3)
            error(" :problem error: SHALLOW WATER equations can only be solved on 1D and 2D grids!")
        end
        inputs[:nvars] = nvars
        println( " # nvars     ", nvars)
        
    elseif (lowercase(inputs[:problem]) == "ns")
        inputs[:problem] = ns()
        
        if (inputs[:nsd] == 1)
            nvars = 3
        elseif(inputs[:nsd] == 2)
            nvars = 4
        elseif(inputs[:nsd] == 3)
            nvars == 5
        end
        inputs[:nvars] = nvars
        println( " # nvars     ", nvars)
        
    elseif (lowercase(inputs[:problem]) == "linearclaw()" ||
            lowercase(inputs[:problem]) == "linclaw" ||
            lowercase(inputs[:problem]) == "lclaw")
        inputs[:problem] = LinearCLaw()
        
        inputs[:nvars] = nvars = 3
        println( " # nvars     ", nvars)
        
    elseif (lowercase(inputs[:problem]) == "advdiff" ||
        lowercase(inputs[:problem]) == "advdif" ||
        lowercase(inputs[:problem]) == "ad" ||
        lowercase(inputs[:problem]) == "adv2d")
        inputs[:problem] = AdvDiff()
        
        inputs[:nvars] = nvars = 1
        println( " # nvars     ", nvars)
    else
        
        inputs[:nvars] = 1 #default
        
        s = """
                jexpresso  user_inputs.jl: problem ", inputs[:problem, " is not coded!
                Chose among:
                    - "AdvDiff"/"AD"/"Adv"
                    - "LinearCLaw"/"LinClaw"
                    - "NS"
                    - "SW"
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
