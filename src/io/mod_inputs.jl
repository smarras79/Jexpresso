using Crayons.Box
using PrettyTables

#macro datatype(str); :($(Symbol(str))); end

function mod_inputs_user_inputs!(parsed_equations, parsed_equations_case_name, equations_dir::String)

    error_flag::Int8 = 0
    
    #
    # Notice: we need `@Base.invokelatest` to call user_inputs() because user_inputs()
    # was defined within this same function via the include(input_dir) above.
    # 
    input_dir = string(equations_dir, "/", parsed_equations, "/", parsed_equations_case_name, "/user_inputs.jl")
    include(input_dir)
    inputs = @Base.invokelatest(user_inputs())

    #Store parsed arguments xxx into inputs[:xxx]
    _parsedToInputs(inputs, parsed_equations, parsed_equations_case_name)
    
    #
    print(GREEN_FG(string(" # Read inputs dict from ", input_dir, " ... \n")))
    pretty_table(inputs; sortkeys=true, border_crayon = crayon"yellow")    
    print(GREEN_FG(string(" # Read inputs dict from ", input_dir, " ... DONE\n")))
    
    #
    # Check that necessary inputs exist in the Dict inside .../IO/user_inputs.jl
    #
    mod_inputs_check(inputs, :nop, Int8(4), "w")  #Polynomial order

    #
    # Plotting parameters:
    #
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

    # Write png to surface using Spline2D interpolation of unstructured data:
    if(!haskey(inputs, :lplot_surf3d))
        inputs[:lplot_surf3d] = false
    end
    if(!haskey(inputs, :smoothing_factor))
        #This is the spline2d smoothing factor. Too small and it may break the spline2d, but it should be as small as possible for precision
        inputs[:smoothing_factor] = 1.0e-1
    end
    
    #
    # END Plotting parameters:
    #
    
    
    #Time:
    if(!haskey(inputs, :ndiagnostics_outputs) && !haskey(inputs, :ndiagnostics_output))
        inputs[:ndiagnostics_outputs] = 2
        inputs[:ndiagnostics_output]  = 2
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
    # https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/
    #
    if(!haskey(inputs, :ode_solver))
        s = """
                        WARNING in user_inputs.jl --> :ode_solver
                        
                            See usable solvers at
                            https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/

                        SSPRK53 will be used by default.
                            """            
            inputs[:ode_solver] = SSPRK54()
        
            @warn s
    end
    
    if(!haskey(inputs, :output_dir))
        inputs[:output_dir] = ""
    end
    if(!haskey(inputs, :loutput_pert))
        inputs[:loutput_pert] = false
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
        inputs[:νx] = Float64(0.0) #default kinematic viscosity
    end
    if(!haskey(inputs, :νy))
        inputs[:νy] = Float64(0.0) #default kinematic viscosity
    end
    if(!haskey(inputs, :νz))
        inputs[:νz] = Float64(0.0) #default kinematic viscosity
    end
    if(!haskey(inputs, :νρ))
        inputs[:νρ] = Float64(0.0) #default kinematic viscosity
    end
    if(!haskey(inputs, :κ))
        inputs[:κ] = Float64(0.0) #default kinematic viscosity
    end

    #
    # Viscous models:
    #
    if(!haskey(inputs, :lvisc))
        inputs[:lvisc] = false
    end
    if(!haskey(inputs, :visc_model))
        inputs[:visc_model] = "av" #Default is artificial viscosity with constant coefficient
    else
        inputs[:visc_model] = lowercase(inputs[:visc_model])
    end
    
    
    #
    # Correct quantities based on a hierarchy of input variables
    #
    # Define default npx,y,z for native grid given
    # values for the user's nelx,y,z
    if(haskey(inputs, :nelx))
        inputs[:npx] = inputs[:nelx] + 1
    else
        inputs[:npx] = UInt8(2)
    end
    if(haskey(inputs, :nely))
        inputs[:npy] = inputs[:nely] + 1
    else
        inputs[:npy] = UInt8(2)
    end
    if(haskey(inputs, :nelz))
        inputs[:npz] = inputs[:nelz] + 1
    else
        inputs[:npz] = UInt8(2)
    end
    
    if (inputs[:nsd] == 1)
        inputs[:npy] = UInt8(1)
        inputs[:npz] = UInt8(1)
    elseif(inputs[:nsd] == 2)
        inputs[:npz] = UInt8(1)
    end

    #Penalty constant for SIPG
    if(!haskey(inputs, :penalty))
        inputs[:penalty] = Float16(0.0) #default kinematic viscosity
    end
    
    
    #------------------------------------------------------------------------
    #To add a new set of governing equations, add a new equations directory
    #to src/equations and call it `ANY_NAME_YOU_WANT` 
    #and add the following lines 
    #
    #elseif (lowercase(equations) == "ANY_NAME_YOU_WANT")
    #inputs[:equations] = ANY_NAME_YOU_WANT()
    #
    #neqs = INTEGER VALUE OF THE NUMBER OF UNKNOWNS for this equations.
    #prinetln( " # neqs     ", neqs)
    #end
    #------------------------------------------------------------------------
    
    #------------------------------------------------------------------------
    # Define neqs based on the equations being solved
    #------------------------------------------------------------------------
    neqs::Int8 = 1
    if (lowercase(parsed_equations) == "compeuler")
        inputs[:equations] = CompEuler()
        inputs[:ldss_laplace] = false
        inputs[:ldss_differentiation] = false
    elseif (lowercase(parsed_equations) == "burgers")
        inputs[:equations] = Burgers()
        inputs[:ldss_laplace] = false
        inputs[:ldss_differentiation] = false
    elseif (lowercase(parsed_equations) == "shallowwater")
        inputs[:equations] = ShallowWater()    
        inputs[:ldss_laplace] = false
        inputs[:ldss_differentiation] = false    
    elseif (lowercase(parsed_equations) == "advdiff" ||
        lowercase(parsed_equations) == "advdif" ||
        lowercase(parsed_equations) == "ad" ||
        lowercase(parsed_equations) == "adv2d")
        inputs[:equations] = AdvDiff()
        inputs[:ldss_laplace] = false
        inputs[:ldss_differentiation] = false
    elseif (lowercase(parsed_equations) == "elliptic" ||
        lowercase(parsed_equations) == "diffusion")
        inputs[:equations] = Elliptic()
        inputs[:ldss_laplace] = true
        inputs[:ldss_differentiation] = false     
    else
        
        #inputs[:neqs] = 1 #default
        
        s = """
                jexpresso  user_inputs.jl: equations ", the inputs[:equations] " that you chose is not coded!
                Chose among:
                         - "CompEuler"
                         - "AdvDiff"/"AD"/"Adv"
              """
        
        @error s
    end

    if(!haskey(inputs, :energy_equation))
        inputs[:energy_equation] = "theta"
        inputs[:δtotal_energy] = 0.0
    else
        if (lowercase(inputs[:equation_set]) == "totalenergy" ||
            lowercase(inputs[:equation_set]) == "totalene"    ||
            lowercase(inputs[:equation_set]) == "totene"      ||
            lowercase(inputs[:equation_set]) == "tene")
            inputs[:δtotal_energy] = 1.0
        else
            #Default
            inputs[:energy_equation] = "theta"
            inputs[:δtotal_energy] = 0.0
        end
    end
    if(!haskey(inputs, :case))
        inputs[:case] = ""
    else
        inputs[:case] = lowercase(inputs[:case])
    end
    if(!haskey(inputs, :lsource))
        inputs[:lsource] = false
    end

    if(!haskey(inputs, :ldss_differentiation))
        inputs[:ldss_differentiation] = false
    end
    if(!haskey(inputs, :ldss_laplace))
        inputs[:ldss_laplace] = false
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


    return inputs
end


function _parsedToInputs(inputs, parsed_equations, parsed_equations_case_name)
    #
    # USER: DO NOT MODIFY inputs[:parsed_equations] and inputs[:parsed_equations_case_name]
    #
    inputs[:parsed_equations]           = parsed_equations
    inputs[:parsed_equations_case_name] = parsed_equations_case_name
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
