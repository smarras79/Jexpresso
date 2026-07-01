using Crayons.Box
using PrettyTables

function mod_inputs_user_inputs!(inputs, rank = 0)

    error_flag::Int8 = 0

    #Store parsed arguments xxx into inputs[:xxx]
    _parsedToInputs(inputs, parsed_equations, parsed_equations_case_name)

    #---------------------------------------------------------------------------
    # Element-Learning pipeline overrides (JEXPRESSO_EL_* environment variables).
    #
    # The automated element-learning pipeline (tools/EL_training/run_element_learning.sh)
    # drives the SAME case through two SEM passes in separate Julia processes:
    #
    #   Phase 1 (sample)     JEXPRESSO_EL_SAMPLE=true   + the 1x1 (single-element) mesh
    #                        → writes input_tensor.csv / output_tensor.csv
    #   Phase 3 (inference)  JEXPRESSO_EL_SAMPLE=false  + the NxN (multi-element) mesh
    #                        → loads the trained :NNfile and writes the solution
    #
    # Reading these knobs from the environment lets the launcher flip lEL_Sample,
    # the mesh, the model file and the output directory per phase WITHOUT editing
    # (and re-committing) the case's user_inputs.jl. Any variable that is unset or
    # empty leaves the user_inputs.jl value untouched, so a plain `run_case` is
    # completely unaffected.
    #---------------------------------------------------------------------------
    _el_env_override!(inputs, rank)

    print_rank(GREEN_FG(string(" # Read inputs dict from ", user_input_file, " ... \n")); msg_rank = rank)
    if rank == 0
        # Wrap long values across multiple lines so nothing is cropped — the
        # default pretty_table truncates wide cells with "…" and hides the
        # tail of long paths / vectors / NamedTuples that users need to see.
        term_cols = try displaysize(stdout)[2] catch; 120 end
        key_w     = 32
        # Leave room for the two outer borders, the column separator, and the
        # padding PrettyTables adds around each cell (≈ 7 chars total).
        val_w     = max(40, term_cols - key_w - 7)
        pretty_table(inputs;
                     sortkeys       = true,
                     border_crayon  = crayon"yellow",
                     linebreaks     = true,
                     autowrap       = true,
                     columns_width  = [key_w, val_w],
                     crop           = :none)
    end
    print_rank(GREEN_FG(string(" # Read inputs dict from ", user_input_file, " ... DONE\n")); msg_rank = rank)
    
    #
    # Check that necessary inputs exist in the Dict inside .../IO/user_inputs.jl
    #
    mod_inputs_check(inputs, :nop, Int8(4), "w")  #Polynomial order
    
    if(!haskey(inputs, :backend))
        inputs[:backend] = CPU()
    end
    
    if (inputs[:backend] != CPU())
        if (inputs[:backend] == CUDABackend())
            global TInt = Int32
            global TFloat = Float32
            global cpu = false
        else
            global TInt = Int32
            global TFloat = Float32
            global cpu = false
        end
    end

    if(!haskey(inputs, :lmanufactured_solution))
       inputs[:lmanufactured_solution] = false
    end

    if(!haskey(inputs, :RT_amr_threshold))
       inputs[:lRT_amr_threshold] = [1.0]
    end

    if(!haskey(inputs, :lRT_problem))
       inputs[:lRT_problem] = false
    end

    if(!haskey(inputs, :lRT_from_data))
       inputs[:lRT_from_data] = false
    end

    if(!haskey(inputs, :RT_data_file))
       inputs[:RT_data_file] = ""
    end

    if(!haskey(inputs, :lcubed_sphere_angular_mesh))
       inputs[:lcubed_sphere_angular_mesh] = false
    end

    if(!haskey(inputs, :rad_HG_g))
      inputs[:rad_HG_g] = 0
    end

    if(!haskey(inputs, :extra_dimensions))
      inputs[:extra_dimensions] = 0
    end
    
    if(!haskey(inputs, :adaptive_extra_meshes))
      inputs[:adaptive_extra_meshes] = false
    end

    if(!haskey(inputs, :extra_dimensions_order))
      inputs[:extra_dimensions_order] = 0
    end

    if(!haskey(inputs, :extra_dimensions_nelemx))
      inputs[:extra_dimensions_nelemx] = 2
    end

    if(!haskey(inputs, :extra_dimensions_nelemy))
      inputs[:extra_dimensions_nelemy] = 2
    end

    if(!haskey(inputs, :extra_dimensions_xmin))
      inputs[:extra_dimensions_xmin] = 0
    end

    if(!haskey(inputs, :extra_dimensions_xmax))
        inputs[:extra_dimensions_xmax] = 2*π
    end

    if(!haskey(inputs, :extra_dimensions_ymin))
      inputs[:extra_dimensions_ymin] = 0
    end

    if(!haskey(inputs, :extra_dimensions_ymax))
        inputs[:extra_dimensions_ymax] = 2*π
    end

    if(!haskey(inputs, :lwall_model))
       inputs[:lwall_model] = false
    end

    # Default to the rank-0-read + MPI.bcast mesh path on every platform.
    # The alternative — `GmshDiscreteModel(parts, file)` in the
    # `lxy_partition=false` branch of `mod_mesh_read_gmsh!` — goes
    # through GridapGmsh's "distributed" constructor, which (depending
    # on the release) parses the .msh file on every rank: nparts × file
    # I/O, nparts × gmsh parses, nparts × peak GMSH memory. On a laptop
    # with a non-trivial mesh that adds minutes to pre-processing
    # before the time-loop even starts.
    #
    # Originally this was macOS-only because the parallel constructor
    # SIGBUSes on Apple Silicon + Open MPI. The serial-read + bcast
    # path has since been the macOS default with no issues, so make it
    # the default on Linux too. Users who need a different partition
    # strategy can still opt out by setting `:lxy_partition => false`
    # in their user_inputs.jl.
    if(!haskey(inputs, :lxy_partition))
        inputs[:lxy_partition] = true
    end

    if(!haskey(inputs, :ifirst_wall_node_index))
         inputs[:ifirst_wall_node_index] = 2 #default is the first LGL point above the surface node along the vertical direction of the surface element
    end
    
    if(!haskey(inputs, :lkep))
       inputs[:lkep] = false
    end
    
    if(!haskey(inputs, :bdy_fluxes))
       inputs[:bdy_fluxes] = false
    end

    if(!haskey(inputs, :bulk_fluxes))
        inputs[:bulk_fluxes] = false
    else
        if inputs[:bulk_fluxes]  == true
            if inputs[:bdy_fluxes]  == false
                inputs[:bdy_fluxes]  = true
            end
        end
    end
    
    if(!haskey(inputs, :LST))
       inputs[:LST] = false
    end

    if(!haskey(inputs, :LST_files))
        inputs[:LST_files] = ("./data_files/LS_heat_forcing.dat","./data_files/LS_rad_cooling.dat","./data_files/LS_vapor_forcing.dat")
    end

    if(!haskey(inputs, :nlay_pg))
       inputs[:nlay_pg] = 10
    end

    if(!haskey(inputs, :nx_pg))
       inputs[:nx_pg] = 10
    end

    if(!haskey(inputs, :ny_pg))
       inputs[:ny_pg] = 10
    end

    if(!haskey(inputs, :ltwo_stream_radiation))
       inputs[:ltwo_stream_radiation] = false
    end

    if(!haskey(inputs, :lphysics_grid))
       inputs[:lphysics_grid] = false
    end

    if(!haskey(inputs, :sounding_file))
       inputs[:sounding_file] = "empty"
    end

    if(!haskey(inputs, :topo_database))
       inputs[:topo_database] = "empty"
    end

    if(!haskey(inputs, :read_topo_latmin))
        inputs[:read_topo_latmin] = -89.99
    end

    if(!haskey(inputs, :read_topo_latmax))
        inputs[:read_topo_latmax] = 89.99
    end
    
    if(!haskey(inputs, :read_topo_lonmin))
        inputs[:read_topo_lonmin] = -179.99
    end

    if(!haskey(inputs, :read_topo_lonmax))
        inputs[:read_topo_lonmax] = 179.99
    end

    if(!haskey(inputs, :read_topo_zone))
        inputs[:read_topo_zone] = 20
    end

    if(!haskey(inputs, :llinsolve))
      inputs[:llinsolve] = false
    end

    if(!haskey(inputs, :lsparse))
        inputs[:lsparse] = true

        if(haskey(inputs, :lelementLearning) &&
            inputs[:lelementLearning] == true)
           # inputs[:lsparse] = false
        end
    else
        if(inputs[:lsparse] == true &&
            haskey(inputs, :lelementLearning) &&
            inputs[:lelementLearning] == true)
            #inputs[:lsparse] = false
        end
    end
    
    if(!haskey(inputs, :NNfile))
      inputs[:NNfile] = nothing
    end


    if(!haskey(inputs, :plot_vlines))
      inputs[:plot_vlines] = "empty"
    end

    if(!haskey(inputs, :plot_hlines))
      inputs[:plot_hlines] = "empty"
    end

    # Colormap for the 2D PNG writer (any ColorSchemes.jl name). The
    # default is cmocean's desaturated diverging "balance", which renders
    # wave fields better than highly saturated maps like viridis.
    if(!haskey(inputs, :plot_colormap))
      inputs[:plot_colormap] = :balance
    end

    # PNG writers: true (default) renders all variables of an output time
    # as ONE plot-matrix figure -- the gksqt window is updated in place and
    # fields-it<n>.png is written. false writes one silent PNG per variable
    # instead (<var>-it<n>.png) and opens no window; see render_plot_matrix
    # in jeplots.jl for why the two modes are mutually exclusive under GR.
    if(!haskey(inputs, :plot_matrix))
      inputs[:plot_matrix] = true
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

    if(!haskey(inputs,:mu_z))
        inputs[:mu_z] = 0.0
    end

    if(!haskey(inputs,:lwarp))
        inputs[:lwarp] = false
    end

    if inputs[:lwarp] == true
        if(!haskey(inputs,:z_transition_start))
            inputs[:z_transition_start] = -1000.0
        end
        if(!haskey(inputs,:z_transition_end))
             inputs[:z_transition_end] = 2200.0
        end
    end
    
    if(!haskey(inputs,:lstretch))
        inputs[:lstretch] = false
    end
    
    if inputs[:lstretch] == true
        if(!haskey(inputs,:stretch_type))
            inputs[:stretch_type] = "powerlaw"
        else
            if inputs[:stretch_type] == "fixed_first"
                if(!haskey(inputs, :first_zelement_size))
                    inputs[:first_zelement_size] = 1.0;
                end
            elseif (inputs[:stretch_type] == "fixed_first_twoblocks_weak" ||
                inputs[:stretch_type] == "fixed_first_twoblocks_strong" ||
                inputs[:stretch_type] == "fixed_first_twoblocks_strong_weak")
                
                if(!haskey(inputs, :first_zelement_size))
                    inputs[:first_zelement_size] = 1.0;
                end
                if(!haskey(inputs, :max_zelement_size_bottom))
                    inputs[:max_zelement_size_bottom] = 1.0;
                end
                if(!haskey(inputs, :zlevel_transition))
                    inputs[:zlevel_transition] = 1000000000.0
                end
                if(!haskey(inputs, :uniform_zelement_size))
                    inputs[:uniform_zelement_size] = 1.0
                end
                if(!haskey(inputs, :max_zelement_size_top))
                    inputs[:max_zelement_size_top] = 1.0;
                end
            end
            
        end
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
        elseif lowercase(inputs[:outformat]) == "netcdf" || lowercase(inputs[:outformat]) == "netcdf"
            inputs[:outformat] = NETCDF()
        end
    end

    # Write png to surface using Spline2D interpolation of unstructured data:
    if(!haskey(inputs, :lplot_surf3d))
        inputs[:lplot_surf3d] = false
    end
    if(!haskey(inputs, :lvolume3d))
        inputs[:lvolume3d] = false
    end
    if(!haskey(inputs, :smoothing_factor))
        #This is the spline2d smoothing factor. Too small and it may break the spline2d, but it should be as small as possible for precision
        inputs[:smoothing_factor] = 1.0e-1
    end
    
    #
    # END Plotting parameters:
    #

    #Restart:
    if (!haskey(inputs, :lrestart))
        inputs[:lrestart] = false
    end
    #
    # Time:
    #
    if(!haskey(inputs, :ndiagnostics_outputs))
        inputs[:ndiagnostics_outputs] = 0
    end
    if(!haskey(inputs, :Δt))
        inputs[:Δt] = 0.1  #Initial time is 0.0 by default
    end
    
    if(!haskey(inputs, :radiation_time_step))
        inputs[:radiation_time_step] = inputs[:Δt]*100
    end

    if(!haskey(inputs, :restart_time))
        inputs[:restart_time] = 0.0
    end

    #mod_inputs_check(inputs, :Δt, Float64(0.1), "w") #Δt --> this will be computed from CFL later on
    if(!haskey(inputs, :tinit))
        inputs[:tinit] = 0.0  #Initial time is 0.0 by default
    end
    if(!haskey(inputs, :tend))
        inputs[:tend] = 0.0  #end time is 0.0 by default
    end

    if( !haskey(inputs, :diagnostics_at_times) )
        inputs[:diagnostics_at_times] = inputs[:tend]
        if (!haskey(inputs, :ndiagnostics_outputs))
            inputs[:ndiagnostics_outputs] = 1 #Force this to none to avoid double output
        end
    else
        inputs[:ndiagnostics_outputs] = 0
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
    # Element learning (lelemLearning)
    #
    if (!haskey(inputs, :Nsamp))
    	inputs[:Nsamp] = 1
    end
    if (!haskey(inputs, :lelementLearning))
        inputs[:lelementLearning] = false
    end

    # Classical FFT (Fourier spectral) Laplace/Poisson solver — selected in
    # problems/drivers.jl when :llinsolve is on. Off by default.
    if (!haskey(inputs, :lfft))
        inputs[:lfft] = false
    end

    # Chebyshev spectral (collocation) Laplace/Poisson solver — the non-periodic
    # counterpart of :lfft. Off by default.
    if (!haskey(inputs, :lcheb))
        inputs[:lcheb] = false
    end

    if !haskey(inputs, :lEL_Sample)
        inputs[:lEL_Sample] = false
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
        inputs[:lwrite_initial] = true
    end

    if (!haskey(inputs, :gmsh_filename_c))
        if haskey(inputs, :gmsh_filename)
            inputs[:gmsh_filename_c] = inputs[:gmsh_filename]
        else
            inputs[:gmsh_filename_c] = "none"
        end
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
        mod_inputs_check(inputs, :gmsh_filename_c, "e")
        
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

    
    if (!haskey(inputs, :lwarmup))
        inputs[:lwarmup] = false
    else
        if !haskey(inputs, :gmsh_filename_c)
            if haskey(inputs, :gmsh_filename)
                inputs[:gmsh_filename_c] = inputs[:gmsh_filename]
            else
                inputs[:gmsh_filename_c] = "none"
            end
        end
    end
    
    #
    # Some physical constants and parameters:
    #
    if(!haskey(inputs, :μ))
        inputs[:μ] = (Float64(0.0)) #default kinematic viscosity
    end

#if(!haskey(inputs, :ivisc_equations))
#        inputs[:ivisc_equations] = [1]
#    end

    #
    # DSGS
    #
    if(!haskey(inputs, :C1))
        inputs[:C1] = 0.0
    end
    if(!haskey(inputs, :C2))
        inputs[:C2] = 0.0
    end
    if(!haskey(inputs, :Pr))
        inputs[:Pr] = 0.7
    end


    #
    # Viscous models:
    #
    if(!haskey(inputs, :lvisc))
        inputs[:lvisc] = false
    end
    if(!haskey(inputs, :visc_model))
        inputs[:visc_model] = AV() #Default is artificial viscosity with constant coefficient
    end

    
    if(!haskey(inputs, :lrichardson))
        inputs[:lrichardson] = false #Default is artificial viscosity with constant coefficient
    end

    #
    # Kinetic Energy or Entropy Preserving
    #
    if(!haskey(inputs, :lkep))
        inputs[:lkep] = false
    end
        
    if inputs[:lkep] == true
        if(!haskey(inputs, :volume_flux))
            inputs[:volume_flux] = ranocha()
        end
    else
        if(!haskey(inputs, :volume_flux))
            inputs[:volume_flux] = nothing
        end
    end

    if(!haskey(inputs, :entropy_variables))
        inputs[:entropy_variables] = false
    end

    #
    # saturation adjustment:
    #
    if(!haskey(inputs, :lsaturation))
        inputs[:lsaturation] = false
    end

    

    #
    # Array of user-defined constant with a user-given meaning. For example, this is used in drivers for the elliptic problems
    #
    if(!haskey(inputs, :rconst))
        inputs[:rconst] = Float64(0.0)
    end
    if(!haskey(inputs, :iconst))
        inputs[:iconst] = Int32(1)
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
    
    if(!haskey(inputs, :lsponge))
        inputs[:lsponge] = false
    end
    if(!haskey(inputs, :zsponge))
        inputs[:zsponge] = 14000.0
    end
    if  inputs[:lsponge] == true
        if(!haskey(inputs, :zsponge))
            inputs[:zsponge] = 14000.0
        end
    end

    if(!haskey(inputs, :lmoist))
        inputs[:lmoist] = false
    end

    if(!haskey(inputs, :lprecip))
        inputs[:lprecip] = false
    end
    
    if(!haskey(inputs, :energy_equation))
        inputs[:energy_equation] = "theta"
        inputs[:δtotal_energy] = 0.0
    end

    if(!haskey(inputs, :lrichardson))
        inputs[:lrichardson] = false
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
    
    if(!haskey(inputs, :loverwrite_output))
        inputs[:loverwrite_output] = false
    end

    if(!haskey(inputs, :SOL_VARS_TYPE))
        inputs[:SOL_VARS_TYPE] = TOTAL() #vs PERT()
    end

    if(!haskey(inputs, :sol_vars_names))
        inputs[:sol_vars_names] = ("q1", "q2", "q3", "q4")
    end
    
    if(!haskey(inputs, :case))
        inputs[:case] = ""
    else
        inputs[:case] = lowercase(inputs[:case])
    end
    if(!haskey(inputs, :lsource))
        inputs[:lsource] = false
    end
    if(!haskey(inputs, :luser_function))
        inputs[:luser_function] = false
    end

    if(!haskey(inputs, :ldss_differentiation))
        inputs[:ldss_differentiation] = false
    end
    if(!haskey(inputs, :ldss_laplace))
        inputs[:ldss_laplace] = false
    end

    # AMR    
    if(!haskey(inputs, :lamr))
        inputs[:lamr] = false
    end

    # LES statistics defaults (used by giga_les TimeIntegrators.jl callbacks).
    if(!haskey(inputs, :statistics_time))
        inputs[:statistics_time] = Float64[]
    end
    if(!haskey(inputs, :statistics_online_start))
        inputs[:statistics_online_start] = Inf
    end
    if(!haskey(inputs, :statistics_online_interval))
        inputs[:statistics_online_interval] = Float32(inputs[:Δt])
    end

    # VTK / AMR restart defaults (used by giga_les TimeIntegrators.jl).
    if(!haskey(inputs, :lrestart_vtk))
        inputs[:lrestart_vtk] = false
    end
    if(!haskey(inputs, :lrestart_amr))
        inputs[:lrestart_amr] = false
    end

    # LES profile/stress var defaults (used by giga_les params_setup.jl).
    if(!haskey(inputs, :lesprofile_vars))
        inputs[:lesprofile_vars] = []
    end
    if(!haskey(inputs, :lesstress_vars))
        inputs[:lesstress_vars] = []
    end

    if(!haskey(inputs, :amr_freq))
        inputs[:amr_freq] = 0
    end

    if(!haskey(inputs, :amr_max_level))
        inputs[:amr_max_level] = 0
    end
     
    if(!haskey(inputs, :ladapt))
        inputs[:ladapt] = false
    end

    if inputs[:lamr] == true
        inputs[:ladapt] = true
    end

    if(!haskey(inputs, :linitial_refine))
        inputs[:linitial_refine] = false
    end
        
    if(!haskey(inputs, :init_refine_lvl))
        inputs[:init_refine_lvl] = 0
    end
        
    if(!haskey(inputs, :amr_max_level))
        inputs[:amr_max_level] = 0
    end

    if(!haskey(inputs, :lpreadapt))
        inputs[:lpreadapt] = false
    end

    if(!haskey(inputs, :preadapt_max_level))
        inputs[:preadapt_max_level] = 0
    end

    if(!haskey(inputs, :amr_start_time))
        inputs[:amr_start_time] = Float32(0.0)
    end

    if(!haskey(inputs, :user_heatflux))
        inputs[:user_heatflux] = 0.0
        inputs[:δhf] = 0.0
    else
        inputs[:δhf] = 1.0
    end

    if inputs[:lpreadapt] == true
        inputs[:ladapt] = true
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


function mod_inputs_check(inputs, key, error_or_warning::String)
    
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


function mod_inputs_check(inputs, key, value, error_or_warning::String)

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

function mod_inputs_print_welcome(rank = 0)
    if rank == 0
        print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))
        print(BLUE_FG(" # Welcome to ", RED_FG("jexpresso\n")))
        print(BLUE_FG(" # A Julia code to solve conservation laws with continuous spectral elements\n"))
        print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))
    end

end

#---------------------------------------------------------------------------
# Element-Learning pipeline: apply JEXPRESSO_EL_* environment overrides.
#
# Used by the automated pipeline (tools/EL_training/run_element_learning.sh) to
# drive one case through its sampling and inference phases without editing
# user_inputs.jl. Each variable is optional; an unset/empty value is ignored.
#
#   JEXPRESSO_EL_SAMPLE      "true"/"false"  → inputs[:lEL_Sample]
#   JEXPRESSO_EL_MESH        path            → inputs[:gmsh_filename]
#   JEXPRESSO_EL_NNFILE      path            → inputs[:NNfile]
#   JEXPRESSO_EL_OUTPUT_DIR    path          → inputs[:output_dir]
#   JEXPRESSO_EL_NSAMP         integer       → inputs[:Nsamp]
#   JEXPRESSO_EL_TAG           string        → inputs[:EL_tensor_tag]
#                                              (tags the input/output CSV names)
#   JEXPRESSO_EL_DIAGNOSTICS    "true"/"false"→ inputs[:lEL_diagnostics]
#   JEXPRESSO_EL_TIMING_SECONDS number        → inputs[:EL_timing_seconds]
#                                               (BenchmarkTools time budget/solve)
#---------------------------------------------------------------------------
function _el_env_override!(inputs, rank = 0)
    _getenv(name) = (v = strip(get(ENV, name, "")); isempty(v) ? nothing : String(v))
    _note(msg)    = rank == 0 && print(GREEN_FG(string(" # [EL pipeline] ", msg, "\n")))

    if (v = _getenv("JEXPRESSO_EL_SAMPLE")) !== nothing
        lv = lowercase(v)
        if lv in ("true", "1", "yes", "on")
            inputs[:lEL_Sample] = true
            _note("JEXPRESSO_EL_SAMPLE → :lEL_Sample = true (sampling phase)")
        elseif lv in ("false", "0", "no", "off")
            inputs[:lEL_Sample] = false
            _note("JEXPRESSO_EL_SAMPLE → :lEL_Sample = false (inference phase)")
        else
            @warn "Ignoring unrecognised JEXPRESSO_EL_SAMPLE value (use true/false)" value=v
        end
    end

    if (v = _getenv("JEXPRESSO_EL_MESH")) !== nothing
        inputs[:gmsh_filename] = v
        _note(string("JEXPRESSO_EL_MESH → :gmsh_filename = ", v))
    end

    if (v = _getenv("JEXPRESSO_EL_NNFILE")) !== nothing
        inputs[:NNfile] = v
        _note(string("JEXPRESSO_EL_NNFILE → :NNfile = ", v))
    end

    if (v = _getenv("JEXPRESSO_EL_OUTPUT_DIR")) !== nothing
        inputs[:output_dir] = v
        _note(string("JEXPRESSO_EL_OUTPUT_DIR → :output_dir = ", v))
    end

    if (v = _getenv("JEXPRESSO_EL_NSAMP")) !== nothing
        n = tryparse(Int, v)
        if n === nothing
            @warn "Ignoring non-integer JEXPRESSO_EL_NSAMP" value=v
        else
            inputs[:Nsamp] = n
            _note(string("JEXPRESSO_EL_NSAMP → :Nsamp = ", n))
        end
    end

    if (v = _getenv("JEXPRESSO_EL_TAG")) !== nothing
        inputs[:EL_tensor_tag] = v
        _note(string("JEXPRESSO_EL_TAG → :EL_tensor_tag = ", v,
                     "  (tensors: input_tensor_", v, ".csv / output_tensor_", v, ".csv)"))
    end

    if (v = _getenv("JEXPRESSO_EL_DIAGNOSTICS")) !== nothing
        lv = lowercase(v)
        if lv in ("true", "1", "yes", "on")
            inputs[:lEL_diagnostics] = true
            _note("JEXPRESSO_EL_DIAGNOSTICS → :lEL_diagnostics = true")
        elseif lv in ("false", "0", "no", "off")
            inputs[:lEL_diagnostics] = false
            _note("JEXPRESSO_EL_DIAGNOSTICS → :lEL_diagnostics = false")
        else
            @warn "Ignoring unrecognised JEXPRESSO_EL_DIAGNOSTICS value (use true/false)" value=v
        end
    end

    if (v = _getenv("JEXPRESSO_EL_TIMING_SECONDS")) !== nothing
        s = tryparse(Float64, v)
        if s === nothing
            @warn "Ignoring non-numeric JEXPRESSO_EL_TIMING_SECONDS" value=v
        else
            inputs[:EL_timing_seconds] = s
            _note(string("JEXPRESSO_EL_TIMING_SECONDS → :EL_timing_seconds = ", s))
        end
    end

    #-----------------------------------------------------------------------
    # Auto-tag the sampled tensors with the case name for EVERY element-learning
    # case, so that input_tensor.csv / output_tensor.csv from different EL cases
    # never overwrite one another (e.g. elementLearning_hole → input_tensor_
    # elementLearning_hole.csv). Applied only when the tag was not already set —
    # by JEXPRESSO_EL_TAG above or by an explicit :EL_tensor_tag in the case's
    # user_inputs.jl — so those explicit choices still win. Non-EL cases keep the
    # plain filenames.
    #-----------------------------------------------------------------------
    if !haskey(inputs, :EL_tensor_tag) &&
       (get(inputs, :lelementLearning, false) == true || get(inputs, :lEL_Sample, false) == true)
        casename = String(get(inputs, :_parsed_case_name, ""))
        if !isempty(casename)
            inputs[:EL_tensor_tag] = casename
            _note(string("element-learning case → :EL_tensor_tag = ", casename,
                         "  (tensors: input_tensor_", casename, ".csv / output_tensor_", casename, ".csv)"))
        end
    end

    return inputs
end
