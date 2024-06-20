using Crayons.Box
using PrettyTables

function mod_inputs_user_inputs!(inputs)

    error_flag::Int8 = 0
    
    #Store parsed arguments xxx into inputs[:xxx]
    _parsedToInputs(inputs, parsed_equations, parsed_equations_case_name)
    
    print(GREEN_FG(string(" # Read inputs dict from ", user_input_file, " ... \n")))
    pretty_table(inputs; sortkeys=true, border_crayon = crayon"yellow")    
    print(GREEN_FG(string(" # Read inputs dict from ", user_input_file, " ... DONE\n")))
    
    #
    # Check that necessary inputs exist in the Dict inside .../IO/user_inputs.jl
    #    
    if(!haskey(inputs, :backend))
        inputs[:backend] = CPU()
    end
    
    if (inputs[:backend] != CPU())
        if (inputs[:backend] == CUDABackend())
            global TInt = Int64
            global TFloat = Float64
            global cpu = false
        else
            global TInt = Int32
            global TFloat = Float32
            global cpu = false
        end
    end

    if(!haskey(inputs, :llinsolve))
      inputs[:llinsolve] = false
    end

    if(!haskey(inputs, :plot_vlines))
      inputs[:plot_vlines] = "empty"
    end

    if(!haskey(inputs, :plot_hlines))
      inputs[:plot_hlines] = "empty"
    end
    
    if(!haskey(inputs, :plot_axis))
      inputs[:plot_axis] = "empty"
    end
   
    if(!haskey(inputs, :plot_overlap))
      inputs[:plot_overlap] = false
    end

    if(!haskey(inputs, :lperiodic_1d))
      inputs[:lperiodic_1d] = false
    end
    
    if(!haskey(inputs, :llaguerre_bc))
      inputs[:llaguerre_bc] = false
    end

    if(!haskey(inputs, :laguerre_tag))
      inputs[:laguerre_tag] = "none"
    end

    if(!haskey(inputs, :lperiodic_laguerre))
      inputs[:lperiodic_laguerre] = false
    end

    if(!haskey(inputs,:llaguerre_1d_right))
      inputs[:llaguerre_1d_right] = false
    end

    if(!haskey(inputs,:llaguerre_1d_left))
      inputs[:llaguerre_1d_left] = false
    end

    if(!haskey(inputs,:laguerre_beta))
      inputs[:laguerre_beta] = 1.0
    end
    
    if(!haskey(inputs,:nop_laguerre))
        inputs[:nop_laguerre] = 18
    end
    
    if(!haskey(inputs,:xfac_laguerre))
        inputs[:xfac_laguerre] = 1.0
    end

    if(!haskey(inputs,:yfac_laguerre))
        inputs[:yfac_laguerre] = 1.0
    end
     
    if(!haskey(inputs,:lfilter))
        inputs[:lfilter] = false
    end

    if(!haskey(inputs,:mu_x))
        inputs[:mu_x] = 0.0
    end

    if(!haskey(inputs,:mu_y))
        inputs[:mu_y] = 0.0
    end

    if(!haskey(inputs,:lwarp))
        inputs[:lwarp] = false
    end

    if(!haskey(inputs,:mount_type))
        inputs[:lagnesi] = "agnesi"
    end

    if(!haskey(inputs,:a_mount))
        inputs[:a_mount] = 10000.0
    end

    if(!haskey(inputs,:h_mount))
        inputs[:h_mount] = 1.0
    end
    
    if(!haskey(inputs,:c_mount))
        inputs[:c_mount] = 0.0
    end

    if(!haskey(inputs,:luser_bc))
        inputs[:luser_bc] = true
    end
    
    if(!haskey(inputs,:xscale))
        inputs[:xscale] = 1.0
    end

    if(!haskey(inputs,:yscale))
        inputs[:yscale] = 1.0
    end
    
    if(!haskey(inputs, :xdisp))
        inputs[:xdisp] = 0.0
    end
    
    if(!haskey(inputs, :ydisp))
        inputs[:ydisp] = 0.0
    end

    if(!haskey(inputs, :filter_type))
        inputs[:filter_type] = "erf"
    end
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
        elseif lowercase(inputs[:outformat]) == "hdf5" || lowercase(inputs[:outformat]) == "h5"
            inputs[:outformat] = HDF5()
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

    #Restart:
    if (haskey(inputs, :lrestart))
        if(inputs[:lrestart] == true && !haskey(inputs, :restart_input_file_path))
            mod_inputs_check(inputs, :restart_input_file_path, "e")
        end
    else
        inputs[:lrestart] = false
    end
            
    #Time:
    if(!haskey(inputs, :ndiagnostics_outputs) && !haskey(inputs, :ndiagnostics_output))
        inputs[:ndiagnostics_outputs] = 2
        inputs[:ndiagnostics_output]  = 2
    end
    mod_inputs_check(inputs, :Δt, Float64(0.1), "w") #Δt --> this will be computed from CFL later on
    if(!haskey(inputs, :tinit))
        inputs[:tinit] = 0.0  #Initial time is 0.0 by default
    end
    if(!haskey(inputs, :tend))
        inputs[:tend] = 0.0  #end time is 0.0 by default
    end
    
    if(!haskey(inputs, :lexact_integration))
        inputs[:lexact_integration] = false #Default integration rule is INEXACT
    end
    if(!haskey(inputs, :llump))
        inputs[:llump] = false #Default no-mass lumping (this is only useful if we use Exact integration)
    end
    
    if(haskey(inputs, :interpolation_nodes))
        
        if(lowercase(inputs[:interpolation_nodes]) == "llg"  ||
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
    if(!haskey(inputs, :ode_adaptive_solver))
        inputs[:ode_adaptive_solver] = false
    end
    if(!haskey(inputs, :output_dir))
        inputs[:output_dir] = "none"
    end
    if(!haskey(inputs, :loutput_pert))
        inputs[:loutput_pert] = false
    end
    if(!haskey(inputs, :lwrite_initial))
        inputs[:lwrite_initial] = false
    end

    #Grid entries:
    if(!haskey(inputs, :lread_gmsh) || inputs[:lread_gmsh] == false)
        mod_inputs_check(inputs, :nsd,  Int8(1), "e")

        if (inputs[:nsd] == 1)
            if (!haskey(inputs, :nelx))
                if (haskey(inputs, :npx))
                    inputs[:nelx] = inputs[:npx] - 1
                else
                    @error("jexpresso: :nelx or :npx is missing in .../IO/user_inputs.jl. Set either one!")
                end
            else
                inputs[:npx] = inputs[:nelx] + 1
            end

        elseif (inputs[:nsd] == 2)
            if (!haskey(inputs, :nelx))
                if (haskey(inputs, :npx))
                    inputs[:nelx] = inputs[:npx] - 1
                else
                    @error("jexpresso -->  :nelx or :npx is missing in .../IO/user_inputs.jl. Set either one!")
                end
            else
                inputs[:npx] = inputs[:nelx] + 1
            end

            if (!haskey(inputs, :nely))
                if (haskey(inputs, :npy))
                    inputs[:nely] = inputs[:npy] - 1
                else
                    @error("jexpresso -->  :nelx or :npx is missing in .../IO/user_inputs.jl. Set either one!")
                end
            else
                inputs[:npy] = inputs[:nely] + 1
            end

            #set all z variables to being equal to y to keep the convention x-y where y is the vertical in 2D
            inputs[:zmin] = inputs[:ymin]
            inputs[:zmax] = inputs[:ymax]
            inputs[:npz]  = inputs[:npy]
            inputs[:nelz]  =inputs[:nely]
            
        end
        mod_inputs_check(inputs, :xmin, "e")
        mod_inputs_check(inputs, :xmax, "e")
        mod_inputs_check(inputs, :ymax, "e")
        mod_inputs_check(inputs, :ymax, "e")
                
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

        if (haskey(inputs, :AD) && inputs[:AD] == FD())
            mod_inputs_check(inputs, :npx, "e")
            mod_inputs_check(inputs, :npy, "e")
        end
        
        s= string("jexpresso: Some undefined (but unnecessary) user inputs 
                                  MAY have been given some default values.
                                  The user may want to assess whether some should be explicitly set in user_inputs.jl")
        
        #@warn s
        
    end #lread_gmsh
    #
    # Some physical constants and parameters:
    #
    if(!haskey(inputs, :μ))
        inputs[:μ] = (Float64(0.0)) #default kinematic viscosity
    end

    if(!haskey(inputs, :ivisc_equations))
        inputs[:ivisc_equations] = (1)
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
    # BC
    #
    if(!haskey(inputs, :luser_bc))
        inputs[:luser_bc] = false
    end
    if(!haskey(inputs, :lneumann))
        inputs[:lneumann] = false
    end

    #
    # Correct quantities based on a hierarchy of input variables
    #
    # Define default npx,y,z for native grid given
    # values for the user's nelx,y,z
    if(haskey(inputs, :nelx))
        inputs[:npx] = inputs[:nelx] + 1
    end
    if(haskey(inputs, :nely))
        inputs[:npy] = inputs[:nely] + 1
    end
    if(haskey(inputs, :nelz))
        inputs[:npz] = inputs[:nelz] + 1
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
    elseif (lowercase(parsed_equations) == "helmholtz" ||
        lowercase(parsed_equations) == "diffusion")
        inputs[:equations] = Helmholtz()
        inputs[:ldss_laplace] = true
        inputs[:ldss_differentiation] = false
    else
        
        #inputs[:neqs] = 1 #default
        
        s = """
                jexpresso  user_inputs.jl: equations ", the inputs[:equations] " that you chose is not coded!
                Chose among:
                         - "CompEuler"
                         - "AdvDiff"
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
    if(!haskey(inputs, :CL))
        # :CL stands for Conservation Law.
        # :CL => CL()  means that we solve dq/dt + \nabla.F(q) = S(q)
        # :CL => NCL() means that we solve dq/dt + u.\nabla(q)= S(q)        
        inputs[:CL] = CL()
    end

    if(!haskey(inputs, :AD))
        inputs[:AD] = ContGal()
    else
        if inputs[:AD] != ContGal() && inputs[:AD] != FD()
            @mystop(" :AD can only be either ContGal() or FD() at the moment.")
        end
    end

    if(!haskey(inputs, :nop))
        if (inputs[:AD] == FD())
            inputs[:nop] = 1
        else
            inputs[:nop] = 4
        end
    else
        if (inputs[:AD] == FD())
            inputs[:nop] = 1
        end
    end
    
    if(!haskey(inputs, :loverwrite_output))
        inputs[:loverwrite_output] = false
    end

    if(!haskey(inputs, :SOL_VARS_TYPE))
        inputs[:SOL_VARS_TYPE] = TOTAL() #vs PERT()
    end

    if(!haskey(inputs, :sol_vars_names))
        inputs[:sol_vars_names] = ("rho", "rho.u", "rho.v", "rho.theta")
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
    if inputs[:μ] != (0.0)
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
              jexpresso: $key is missing in problems/equations/PROBLEM_NAME/PROBLEM_CASE_NAME/user_inputs.jl
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
    print(BLUE_FG(" # A Julia code to solve conservation laws with continuous spectral elements\n"))
    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))

end
