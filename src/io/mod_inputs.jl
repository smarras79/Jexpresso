using ArgParse
using Crayons.Box
using PrettyTables
using Revise

export mod_inputs_user_inputs
export mod_inputs_print_welcome

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--opt1"
        help = "an option with an argument"
        "--opt2", "-o"
        help = "another option with an argument"
        arg_type = Int
        default = 0
        "--flag1"
        help = "an option without argument, i.e. a flag"
        action = :store_true
        "arg1"
        help = "a positional argument"
        required = true
    end

    return parse_args(s)
end


function mod_inputs_user_inputs!(problem_name, problem_dir::String)

    error_flag::Int8 = 0
    
    #
    # Notice: we need `@Base.invokelatest` to call user_inputs() because user_inputs()
    # was definied within this same function via the include(input_dir) above.
    # 
    input_dir = string(problem_dir, "/", problem_name, "/user_inputs.jl")
    include(input_dir)
    inputs = @Base.invokelatest(user_inputs())
    
    #
    print(GREEN_FG(string(" # Read inputs dict from ", input_dir, " ... \n")))
    pretty_table(inputs; sortkeys=true, border_crayon = crayon"yellow")    
    print(GREEN_FG(string(" # Read inputs dict from ", input_dir, " ... DONE\n")))
    
    #
    # Check that necessary inputs exist in the Dict inside .../IO/user_inputs.jl
    #
    mod_inputs_check(inputs, :nop, Int8(4), "w")  #Polynomial order

    
    if(!haskey(inputs, :outformat))
        inputs[:outformat] = ASCII()
    else
        if lowercase(inputs[:outformat]) == "png"
            inputs[:outformat] = PNG()
        elseif lowercase(inputs[:outformat]) == "ascii"
            inputs[:outformat] = ASCII()
        elseif lowercase(inputs[:outformat]) == "vtk"
            inputs[:outformat] = VTK()
        end
    end
    
    #Time:
    if(!haskey(inputs, :ndiagnostics_outputs))
        inputs[:ndiagnostics_outputs] = 2
    end
    mod_inputs_check(inputs, :Δt, Float64(1.0), "w") #Δt --> this will be computed from CFL later on
    if(!haskey(inputs, :tinit))
        inputs[:tinit] = 0.0  #Initial time is 0.0 by default
    end
     if(!haskey(inputs, :tend))
        inputs[:tend] = 0.0  #end time is 0.0 by default
    end
    
    if(!haskey(inputs, :lexact_integration))
        inputs[:lexact_integration] = false #Default integration rule is INEXACT
    end

    if(haskey(inputs, :interpolation_nodes))
        
        if(lowercase(inputs[:interpolation_nodes]) == "llg" ||
           lowercase(inputs[:interpolation_nodes]) == "gll" ||
           lowercase(inputs[:interpolation_nodes]) == "lgl")
            inputs[:interpolation_nodes] = LGL()

        elseif(lowercase(inputs[:interpolation_nodes]) == "lg" ||
               lowercase(inputs[:interpolation_nodes]) == "gl")
            inputs[:interpolation_nodes] = LG()
            
        elseif(lowercase(inputs[:interpolation_nodes]) == "cg" ||
               lowercase(inputs[:interpolation_nodes]) == "gc")
            inputs[:interpolation_nodes] = CG()
            
        elseif(lowercase(inputs[:interpolation_nodes]) == "cgl" ||
               lowercase(inputs[:interpolation_nodes]) == "gcl")
            inputs[:interpolation_nodes] = CGL()
        else
            s = """
                ERROR in user_inputs.jl --> :interpolation_nodes
                
                    Chose among:
                     - "lgl"
                     - "lg"
                     - "cg"
                     - "cgl"
              """
    
            error(s)
        end
    else
        #default are LGL
        inputs[:interpolation_nodes] = LGL()
    end

    if(haskey(inputs, :quadrature_nodes))
        
        if(lowercase(inputs[:quadrature_nodes]) == "llg" ||
           lowercase(inputs[:quadrature_nodes]) == "gll" ||
           lowercase(inputs[:quadrature_nodes]) == "lgl")
            inputs[:quadrature_nodes] = LGL()

        elseif(lowercase(inputs[:quadrature_nodes]) == "lg" ||
               lowercase(inputs[:quadrature_nodes]) == "gl")
            inputs[:quadrature_nodes] = LG()
            
        elseif(lowercase(inputs[:quadrature_nodes]) == "cg" ||
               lowercase(inputs[:quadrature_nodes]) == "gc")
            inputs[:quadrature_nodes] = CG()
            
        elseif(lowercase(inputs[:quadrature_nodes]) == "cgl" ||
               lowercase(inputs[:quadrature_nodes]) == "gcl")
            inputs[:quadrature_nodes] = CGL()
        else
            s = """
                ERROR in user_inputs.jl --> :quadrature_nodes
                
                    Chose among:
                     - "lgl"
                     - "lg"
                     - "cg"
                     - "cgl"
              """
            
            error(s)            
        end
    else
        #default are LGL
        inputs[:quadrature_nodes] = LGL()
    end
    
    #
    # DifferentialEquations.jl is used to solved the ODEs resulting from the method-of-lines
    #
    if(haskey(inputs, :ode_solver))
        if(uppercase(inputs[:ode_solver]) == "RK4")
            inputs[:ode_solver] = RK4()
        elseif(uppercase(inputs[:ode_solver]) == "SSPRK22")
            inputs[:ode_solver] = SSPRK22()
        elseif(uppercase(inputs[:ode_solver]) == "SSPRK33")
            inputs[:ode_solver] = SSPRK33()
        elseif(uppercase(inputs[:ode_solver]) == "SSPRK53")
            inputs[:ode_solver] = SSPRK53()
        elseif(uppercase(inputs[:ode_solver]) == "SSPRK54")
            inputs[:ode_solver] =SSPRK54()
        elseif(uppercase(inputs[:ode_solver]) == "SSPRK63")
            inputs[:ode_solver] = SSPRK63()
        elseif(uppercase(inputs[:ode_solver]) == "SSPRK73")
            inputs[:ode_solver] = SSPRK73()
        elseif(uppercase(inputs[:ode_solver]) == "SSPRK104")
            inputs[:ode_solver] = SSPRK104()
        elseif(uppercase(inputs[:ode_solver]) == "CARPENTERKENNEDY2N54")
            inputs[:ode_solver] = CarpenterKennedy2N54()
        elseif(uppercase(inputs[:ode_solver]) == "BICGSTAB" ||
               uppercase(inputs[:ode_solver]) == "BICGSTABLE" ||
               uppercase(inputs[:ode_solver]) == "IterativeSolversJL_BICGSTAB") 
            inputs[:ode_solver] = IterativeSolversJL_BICGSTAB()
        elseif(uppercase(inputs[:ode_solver]) == "GMRES"|| uppercase(inputs[:ode_solver]) == "IterativeSolversJL_GMRES")
            inputs[:ode_solver] = IterativeSolversJL_GMRES()
        else
            s = """
                    WARNING in user_inputs.jl --> :ode_solver
                    
                        See usable solvers at
                        https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/

                    SSPRK53 will be used by default.
                        """            
            inputs[:ode_solver] = SSPRK53()

            @warn s
        end
    end

if(!haskey(inputs, :output_dir))
    inputs[:output_dir] = ""
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

        s= string("jexpresso: Some undefined (but unnecessary) user inputs 
                              MAY have been given some default values.
                              User needs not to worry about them.")
        
        #@warn s
        
    end #lread_gmsh
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

    #Penalty constant for SIPG
    if(!haskey(inputs, :penalty))
        inputs[:penalty] = Float16(0.0) #default kinematic viscosity
    end
    
    
    #------------------------------------------------------------------------
    #To add a new set of governing equations, add a new problem director
    #to src/problems and call it `ANY_NAME_YOU_WANT` 
    #and add the following lines 
    #
    #elseif (lowercase(problem_name) == "ANY_NAME_YOU_WANT")
    #inputs[:problem] = ANY_NAME_YOU_WANT()
    #
    #neqns = INTEGER VALUE OF THE NUMBER OF UNKNOWNS for this problem.
    #prinetln( " # neqns     ", neqns)
    #end
    #------------------------------------------------------------------------
    
    #------------------------------------------------------------------------
# Define neqns based on the problem being solved
#------------------------------------------------------------------------
neqns::Int8 = 1
if (lowercase(problem_name) == "burgers")
    inputs[:problem] = Burgers()
    
    if(inputs[:nsd] == 1)
        neqns = 1
    elseif (inputs[:nsd] == 2)
        neqns = 2
    end
    inputs[:neqns] = neqns
    println( " # Number of equations ", neqns)
    
elseif (lowercase(problem_name) == "sw")
    inputs[:problem] = sw()
    
    if (inputs[:nsd] == 1)
        neqns = 2
    elseif(inputs[:nsd] == 2)
        neqns = 3
    elseif(inputs[:nsd] == 3)
        error(" :problem error: SHALLOW WATER equations can only be solved on 1D and 2D grids!")
    end
    inputs[:neqns] = neqns
    println( " # Number of equations ", neqns)
    
elseif (lowercase(problem_name) == "linearclaw" ||
        lowercase(problem_name) == "linclaw" ||
        lowercase(problem_name) == "lclaw")
    inputs[:problem] = LinearCLaw()
    
    inputs[:neqns] = neqns = 3
    println( " # neqns     ", neqns)
    
elseif (lowercase(problem_name) == "advdiff" ||
        lowercase(problem_name) == "advdif" ||
        lowercase(problem_name) == "ad" ||
        lowercase(problem_name) == "adv2d")
    inputs[:problem] = AdvDiff()
    
    inputs[:neqns] = neqns = 1
    println( " # neqns     ", neqns)
    
elseif (lowercase(problem_name) == "elliptic" ||
        lowercase(problem_name) == "diffusion")
    inputs[:problem] = Elliptic()
    
    inputs[:neqns] = neqns = 1
    println( " # neqns     ", neqns)
    
elseif (lowercase(problem_name) == "helmholtz")
    inputs[:problem] = Helmholtz()
    
    inputs[:neqns] = neqns = 1
    println( " # neqns     ", neqns)
    
else
    
    inputs[:neqns] = 1 #default
    
    s = """
            jexpresso  user_inputs.jl: problem ", inputs[:problem, " is not coded!
            Chose among:
                     - "AdvDiff"/"AD"/"Adv"
                     - "LinearCLaw"/"LinClaw"
                     - "Burgers"
                     - "SW"
          """
    
    @error s
end

#------------------------------------------------------------------------
# The following quantities stored in the inputs[] dictionary are only
# auxiliary and are NEVER to be defined by the user
#------------------------------------------------------------------------
if ((inputs[:νx] != 0.0) || (inputs[:νy] != 0.0) || (inputs[:νz] != 0.0))
    inputs[:δvisc] = 1.0
else
    inputs[:δvisc] = 0.0
end


return inputs, neqns
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
