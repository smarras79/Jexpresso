using Crayons.Box
using PrettyTables

function mod_inputs_user_inputs!(inputs, rank = 0)

    error_flag::Int8 = 0
    
    #Store parsed arguments xxx into inputs[:xxx]
    _parsedToInputs(inputs, parsed_equations, parsed_equations_case_name)
    
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
            global TInt   = Int32
            global TFloat = Float32
            global cpu    = false
        else
            global TInt   = Int32
            global TFloat = Float32
            global cpu    = false
        end
    end


    if(!haskey(inputs, :RT_atmos_coupling))
       inputs[:RT_atmos_coupling] = false
    end

    if (inputs[:RT_atmos_coupling])
        inputs[:RT_radiative_heating] = true
    end


    if(!haskey(inputs, :RT_radiative_heating))
       inputs[:RT_radiative_heating] = false
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
    if (inputs[:lRT_from_data])
       if (!(haskey(inputs, :RT_shortwave)))
            inputs[:RT_longwave] = true
            inputs[:RT_shortwave] = false
       elseif inputs[:RT_shortwave]
            inputs[:RT_longwave] = false
       else
            inputs[:RT_longwave] = true
       end

    end

    if (!(haskey(inputs, :RT_shortwave)))
        inputs[:RT_shortwave] = false
    end
    
    if (!(haskey(inputs, :RT_longwave)))
        inputs[:RT_longwave] = false
    end

    if(!haskey(inputs, :RT_S0_flux))
        inputs[:RT_S0_flux] = 1361.0
    end
    if(!haskey(inputs, :RT_μ0))
        inputs[:RT_μ0] = 0.5
    end
    if(!haskey(inputs, :RT_ϕ0))
        inputs[:RT_ϕ0] = 3*π/4
    end
    if(!haskey(inputs, :RT_δ_beam))
        inputs[:RT_δ_beam] = 0.05
    end
    if(!haskey(inputs, :RT_ϵ_surface))
        inputs[:RT_ϵ_surface] = 0.97
    end
    if(!haskey(inputs, :RT_T_space))
        inputs[:RT_T_space] = 0.0
    end

    if(!haskey(inputs, :RT_data_file))
       inputs[:RT_data_file] = ""
    end

    if(!haskey(inputs, :lcubed_sphere_angular_mesh))
       inputs[:lcubed_sphere_angular_mesh] = false
    end

    if(!haskey(inputs, :rad_HG_g))
      inputs[:rad_HG_g] = 0.0
    end

    if(!haskey(inputs, :extra_dimensions))
      inputs[:extra_dimensions] = 0
    end
    
    if(!haskey(inputs, :adaptive_extra_meshes))
      inputs[:adaptive_extra_meshes] = false
    end

    if(!haskey(inputs, :RT_precond))
      inputs[:RT_precond] = :global_lu
    end

    if(!haskey(inputs, :RT_gmres_tol))
      inputs[:RT_gmres_tol] = 1e-4
    end

    if(!haskey(inputs, :RT_gmres_restart))
      inputs[:RT_gmres_restart] = 100
    end

    if(!haskey(inputs, :RT_asm_ilu_tau))
      inputs[:RT_asm_ilu_tau] = 0.1
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
    # HEVI needs the columnar partition; a fully explicit run does not.
    #
    # The vertically-implicit half of HEVI solves one banded system per
    # VERTICAL COLUMN, so it wants whole columns on a rank. The columnar
    # partition gives exactly that -- measured on rtb_hevi, 0 of 180 local
    # column-instances split at four ranks -- and it also happens to be the
    # partition on which the assembly is correct: every node's ghost-inflation
    # factor is then the same in the assembled RHS as in the mass matrix, so
    # the two cancel.
    #
    # Under the p4est space-filling-curve partition they do not always cancel.
    # Measured on rtb with three ranks: the lumped mass comes out 2.0x its
    # serial value at 1155 nodes per rank and exactly 1.5x at another 330,
    # while the assembled stiffness at those same 330 comes out 1.4965x rather
    # than 1.5x, because a derivative operator's two element contributions are
    # unequal where the two element masses are equal. M^-1 K then stops being
    # skew, which a split scheme cannot survive: the HEVI vertical operator
    # measured max(Re)/max|Im| = +5.07e-02 at four ranks, a +3.5 1/s growth
    # rate, and rtb_hevi died at t ~ 12.
    #
    # So HEVI defaults to the columnar partition and an explicit run keeps
    # p4est. A deck can still say either explicitly; if it asks for HEVI on
    # the SFC partition, build_hevi refuses rather than diverging quietly
    # twelve seconds in.
    #
    # IMEX3D needs it for both halves of the same argument: its preconditioner
    # is the very same per-column banded solve, and the skewness argument above
    # applies to the 3D acoustic operator exactly as it does to the vertical
    # one -- more so, since nothing is left explicit to absorb a growing mode.
    if(!haskey(inputs, :lxy_partition))
        inputs[:lxy_partition] = hevi_enabled(inputs) || imex3d_enabled(inputs)
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
        inputs[:outformat] = NONE()
    else
        if lowercase(inputs[:outformat]) == "png"
            inputs[:outformat] = PNG()
        elseif lowercase(inputs[:outformat]) == "vtk"
            inputs[:outformat] = VTK()
        elseif lowercase(inputs[:outformat]) == "hdf5" || lowercase(inputs[:outformat]) == "h5"
            inputs[:outformat] = HDF5()
        elseif lowercase(inputs[:outformat]) == "netcdf" || lowercase(inputs[:outformat]) == "netcdf"
            inputs[:outformat] = NETCDF()
        else
            inputs[:outformat] = NONE()
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
    
    if (!haskey(inputs, :lrestart_vtk))
        inputs[:lrestart_vtk] = false
    end


    if (!haskey(inputs, :lrestart_amr))
        inputs[:lrestart_amr] = false
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



    #---------------------------------------------------------------------------
    #LES statistics
    #---------------------------------------------------------------------------
    if(!haskey(inputs, :statistics_time))
        inputs[:statistics_time] = Float64[]
    end

    if(!haskey(inputs, :statistics_online_start))
        inputs[:statistics_online_start] = Inf
    end

    if(!haskey(inputs, :statistics_online_interval))
        inputs[:statistics_online_interval] = Float32(inputs[:Δt])
    end

    if(!haskey(inputs, :lesprofile_vars))
        inputs[:lesprofile_vars] = []
    end

    if(!haskey(inputs, :lesstress_vars))
        inputs[:lesstress_vars] = []
    end


    #---------------------------------------------------------------------------
    #END LES statistics
    #---------------------------------------------------------------------------

    
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

    #
    # Grid-only runs and the spherical-shell (2D manifold in 3D) grid.
    #
    #   :lspherical_shell  solve on a CLOSED quadrilateral shell: the manifold
    #                      metrics, RHS and time loop replace the flat ones.
    #                      The GRID itself is read by the ordinary gmsh path in
    #                      mesh.jl, which detects a 2D manifold embedded in 3D
    #                      from the model and keeps z (see `lmanifold` there).
    #   :lgrid_only        build the grid, dump it to VTK, and STOP — no
    #                      initial condition, no time integration. This is the
    #                      switch a user flips while the equations for a new
    #                      geometry are still being written.
    #   :linit_only        one step further: build the grid AND the initial
    #                      condition, write both, and STOP before the time
    #                      integration. :lgrid_only wins if both are set.
    #
    if(!haskey(inputs, :lspherical_shell))
        inputs[:lspherical_shell] = false
    end
    if(!haskey(inputs, :lgrid_only))
        inputs[:lgrid_only] = false
    end
    if(!haskey(inputs, :linit_only))
        inputs[:linit_only] = false
    end
    if(!haskey(inputs, :lcheck_grid))
        inputs[:lcheck_grid] = true
    end
    if(!haskey(inputs, :lstop_on_bad_grid))
        inputs[:lstop_on_bad_grid] = true
    end
    # NOTE :lmerge_coincident_nodes / :node_merge_tol are gone with the bespoke
    # shell reader. A watertight gmsh grid carries one node per seam location
    # already, and Gridap's topology is built from the node ids in the file, so
    # there is nothing to merge. A grid that really does duplicate its seam
    # nodes is a broken grid: fix it in gmsh (share the curves between panels)
    # rather than stitching it back together at read time.
    if(!haskey(inputs, :lproject_to_sphere))
        inputs[:lproject_to_sphere] = true
    end
    #
    #   :cubed_sphere_map         move the shell's nodes onto a different
    #                             cube-face → sphere map after the grid is read.
    #                             Connectivity, panel decomposition and panel
    #                             boundaries are untouched — only where the
    #                             nodes sit inside each panel changes. See
    #                             src/kernel/mesh/cubed_sphere_maps.jl.
    #
    #     :none         (default) leave the grid exactly as the .msh has it.
    #     :gnomonic     equidistant central projection, Sadourny (1972).
    #                   (:equidistant is an accepted alias.)
    #     :equiangular  central projection with the face coordinate measured as
    #                   an angle, Ronchi, Iacono & Paolucci (1996). The most
    #                   HOMOGENEOUS of the three: largest minimum grid distance,
    #                   so the largest explicit time step.
    #     :conformal    Rančić, Purser & Mesinger (1996). ORTHOGONAL AND ISOTROPIC
    #                   everywhere, corner included — and REFUSED on any grid with
    #                   a node on a cube corner, which is every structured panel
    #                   grid, because its Jacobian is zero there. Three panels meet
    #                   at a corner and each opens 120°, so a map that keeps the
    #                   square's 90° can only do it by collapsing its derivative;
    #                   there is no regularisation of it a nodal scheme can use
    #                   (THE 120° CORNER in src/kernel/mesh/cubed_sphere_maps.jl
    #                   has the argument and the measurements). Use :equiangular.
    #                   (:conformal_exact is a deprecated alias.)
    #
    # The grid the remap starts FROM is MEASURED by detect_cubed_sphere_map, not
    # assumed, so there is deliberately no input for it: a "source map" switch is
    # a claim about the .msh that nothing can check, and getting it wrong silently
    # produces a grid that is neither map. (It used to be assumed gnomonic, and
    # that was wrong — gmsh spaces `Transfinite Line` points at equal ANGLE along
    # a `Circle` arc, so cubed_sphere.geo emits the EQUIANGULAR grid.)
    #
    if(!haskey(inputs, :cubed_sphere_map))
        inputs[:cubed_sphere_map] = :none
    end
    let _m = inputs[:cubed_sphere_map]
        (_m === :none || _m in CUBED_SPHERE_MAPS) ||
            error(string(" # ERROR mod_inputs.jl: :cubed_sphere_map => ", _m,
                         " is not recognised. Use :none or one of ", CUBED_SPHERE_MAPS, "."))
    end
    #
    #   :sphere_metrics    the 2D-manifold metric terms of build_sphere_metrics.
    #                      Kopriva's curl-invariant form (J. Sci. Comput. 26(3),
    #                      301, 2006, Eq. 15; Sec. 3.2.3 of Kelly, Alves,
    #                      Eckermann et al., JCP 552, 2026, 114683) DEGENERATES
    #                      on a surface: it fixes the in-surface metric terms
    #                      only up to the direction v in which the surface is
    #                      extended off itself, and both choices below are
    #                      curl-invariant. See the header of sphere_metrics.jl.
    #
    #     :cross_product   (default) v = n̂, the discrete surface normal — which
    #                      is what the textbook cross-product formulas give.
    #                      Uniquely, its strong-form divergence annihilates a
    #                      rigid rotation exactly.
    #     :radial          v = x̂, the exact radial. ∇ₛ is then tangent to the
    #                      TRUE sphere rather than to the interpolant of it, at
    #                      the cost of the rigid-rotation property.
    #                      (:curl_invariant is an accepted alias.)
    #
    if(!haskey(inputs, :sphere_metrics))
        inputs[:sphere_metrics] = :cross_product
    end


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
    # Marras-Nazarov DynSGS (visc_model = DSGS_MHD()) parameters.
    #   :dsgs_C1    coefficient of the residual viscosity  C1·Δ²·‖R‖/‖q−⟨q⟩‖
    #   :dsgs_C2    coefficient of the wave-speed cap      C2·Δ·(|v|+c_f)
    #   :dsgs_gamma ratio of specific heats used by the MHD EOS and the fast
    #               magnetosonic speed (5/3 for the monatomic plasma cases;
    #               deliberately NOT PhysConst.γ, which is air's 1.4)
    #   :dsgs_Prt   turbulent Prandtl number for the energy slot
    #
    if(!haskey(inputs, :dsgs_C1))
        inputs[:dsgs_C1] = 1.0
    end
    if(!haskey(inputs, :dsgs_C2))
        inputs[:dsgs_C2] = 0.5
    end
    if(!haskey(inputs, :dsgs_gamma))
        inputs[:dsgs_gamma] = 5.0/3.0
    end
    if(!haskey(inputs, :dsgs_Prt))
        inputs[:dsgs_Prt] = 0.7
    end
    # Scope of the DynSGS normalising scales ⟨q⟩ and ‖q−⟨q⟩‖.
    #
    #   false (default) : rank-local. No communication at all.
    #   true            : the domain norms of Marras eq. (9) / Nazarov &
    #                     Hoffman eq. (3.5). Costs 2-3 MPI Allreduce per RHS
    #                     call — 10-15 per step under a five-stage RK — on
    #                     every rank's critical path.
    #
    # These two quantities only set the SCALE the element residual is measured
    # against, and a partition of a connected domain resolves that scale as
    # well as the whole domain does, so the default costs nothing and changes
    # the solution only at round-off level. Set it true when μ has to be
    # reproducible across rank counts (bit-for-bit regression tests), or when
    # a rank's subdomain genuinely cannot see the solution's scale. Serial
    # runs are unaffected either way. See kernel/physics/SGS.jl
    # (_dsgs_norm_scope) and ENVIRONMENT_VARIABLES.md.
    if(!haskey(inputs, :ldsgs_global_norms))
        inputs[:ldsgs_global_norms] = false
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

    
    #
    # Dynamic SGS closures (SMAG/VREM). This block is the single source of
    # truth for these keys: params_setup.jl reads them, it does not re-default
    # them.
    #
    if(!haskey(inputs, :lrichardson))
        # Buoyancy (Richardson) correction on the eddy viscosity. Off by
        # default to keep existing decks bit-identical; any stratified case —
        # and every convective boundary layer — wants it ON, otherwise the full
        # eddy diffusivity acts across a capping inversion and smears it.
        inputs[:lrichardson] = false
    end

    if(!haskey(inputs, :C_s))
        # Smagorinsky constant. 0.21 is the historical Jexpresso value; ABL LES
        # normally runs 0.13-0.18 (Lilly 0.17, Deardorff/Moeng ~0.13 with
        # shear). Since nu_t scales with C_s^2, the difference is a factor 1.4-2.6.
        inputs[:C_s] = PhysicalConst{Float64}().C_s
    end

    if(!haskey(inputs, :lwall_damping))
        # Near-wall limit l = min(C_s*Delta, kappa*z) on the SGS mixing length
        # (SMAG only). Off by default: it needs the distance to the wall, which
        # for a warped/terrain-following mesh is not the height above the
        # domain floor.
        inputs[:lwall_damping] = false
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

    # HDF5 solution restart (:lrestart/:restart_time) does not know about the
    # p4est forest and cannot reconstruct the adapted mesh it was written
    # against, so it's incompatible with AMR. Force it off rather than
    # silently corrupting/erroring on a shape mismatch; use :lrestart_amr
    # (VTK + p4est forest restart) for AMR restarts instead.
    if inputs[:lamr] == true
        inputs[:lrestart]     = false
        inputs[:restart_time] = 0.0
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

    if(!haskey(inputs, :amr_start_time))
        inputs[:amr_start_time] = Float32(0.0)
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

function build_tspan(inputs, TFloat)
    if get(inputs, :lamr, false) == true
        amr_freq = inputs[:amr_freq]
        Δt_amr   = amr_freq * inputs[:Δt]
        [TFloat(inputs[:tinit]), TFloat(inputs[:tinit] + inputs[:amr_start_time] + Δt_amr)]
    else
        [TFloat(inputs[:tinit]), TFloat(inputs[:tend])]
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
