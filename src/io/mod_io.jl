using Crayons.Box
using PrettyTables

using ArgParse
function mod_io_parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "eqs"
        help = "Directoy that contains some user-defined cases"
        default = "CompEuler"
        required = false

        "eqs_case"
        help = "case name in equations directory"
        default = "wave1d"
        required = false

        "CI_MODE"
        help = "CI_MODE: true or false"
        default = "false"
        required = false
    end

    return parse_args(s)
end


function mod_io_parse_args()
    
    parsed_args = mod_io_parse_commandline()
    
    # These must be global because other functions access them by name (legacy design)
    global parsed_equations           = string(parsed_args["eqs"])
    global parsed_equations_case_name = string(parsed_args["eqs_case"])
    global parsed_CI_mode             = string(parsed_args["CI_MODE"])
    global driver_file                = string(joinpath( "..", "problems", "drivers.jl"))
    
    if parsed_CI_mode == "true"
        case_name_dir = string(joinpath("test", "CI-runs", parsed_equations, parsed_equations_case_name))
    else
        case_name_dir = string(joinpath("problems", parsed_equations, parsed_equations_case_name))
    end
    
    global user_input_file      = string(joinpath("..", case_name_dir, "user_inputs.jl"))
    global user_flux_file       = string(joinpath("..", case_name_dir, "user_flux.jl"))
    global user_source_file     = string(joinpath("..", case_name_dir, "user_source.jl"))
    global user_bc_file         = string(joinpath("..", case_name_dir, "user_bc.jl"))
    global user_initialize_file = string(joinpath("..", case_name_dir, "initialize.jl"))
    global user_primitives_file = string(joinpath("..", case_name_dir, "user_primitives.jl"))
    
end

# ---------------------------------------------------------------------------
# Apply all simple scalar defaults in one pass (get! = one hash lookup vs two
# for the if(!haskey)…end pattern, and a much smaller compilation unit for
# the Julia JIT).
# ---------------------------------------------------------------------------
function _apply_defaults!(inputs)
    get!(inputs, :backend,                CPU())
    get!(inputs, :lwall_model,            false)
    get!(inputs, :ifirst_wall_node_index, 2)
    get!(inputs, :lkep,                   false)
    get!(inputs, :bdy_fluxes,            false)
    get!(inputs, :bulk_fluxes,           false)
    get!(inputs, :LST,                   false)
    get!(inputs, :LST_files,             ("./data_files/LS_heat_forcing.dat",
                                          "./data_files/LS_rad_cooling.dat",
                                          "./data_files/LS_vapor_forcing.dat"))
    get!(inputs, :nlay_pg,               10)
    get!(inputs, :nx_pg,                 10)
    get!(inputs, :ny_pg,                 10)
    get!(inputs, :ltwo_stream_radiation, false)
    get!(inputs, :lphysics_grid,         false)
    get!(inputs, :sounding_file,         "empty")
    get!(inputs, :topo_database,         "empty")
    get!(inputs, :read_topo_latmin,      -89.99)
    get!(inputs, :read_topo_latmax,       89.99)
    get!(inputs, :read_topo_lonmin,      -179.99)
    get!(inputs, :read_topo_lonmax,       179.99)
    get!(inputs, :read_topo_zone,         20)
    get!(inputs, :llinsolve,             false)
    get!(inputs, :lsparse,               true)
    get!(inputs, :plot_vlines,           "empty")
    get!(inputs, :plot_hlines,           "empty")
    get!(inputs, :plot_axis,             "empty")
    get!(inputs, :plot_overlap,          false)
    get!(inputs, :lperiodic_1d,          false)
    get!(inputs, :llaguerre_bc,          false)
    get!(inputs, :laguerre_tag,          "none")
    get!(inputs, :lperiodic_laguerre,    false)
    get!(inputs, :llaguerre_1d_right,    false)
    get!(inputs, :llaguerre_1d_left,     false)
    get!(inputs, :laguerre_beta,         1.0)
    get!(inputs, :nop_laguerre,          18)
    get!(inputs, :xfac_laguerre,         1.0)
    get!(inputs, :yfac_laguerre,         1.0)
    get!(inputs, :lfilter,               false)
    get!(inputs, :mu_x,                  0.0)
    get!(inputs, :mu_y,                  0.0)
    get!(inputs, :mu_z,                  0.0)
    get!(inputs, :lwarp,                 false)
    get!(inputs, :z_transition_start,    -1000.0)   # meaningful only when lwarp=true
    get!(inputs, :z_transition_end,       2200.0)   # meaningful only when lwarp=true
    get!(inputs, :lwith_alya,            false)
    get!(inputs, :lstretch,              false)
    get!(inputs, :lagnesi,               "agnesi")
    get!(inputs, :a_mount,               10000.0)
    get!(inputs, :h_mount,               1.0)
    get!(inputs, :c_mount,               0.0)
    get!(inputs, :luser_bc,              true)      # first-occurrence default wins
    get!(inputs, :xscale,                1.0)
    get!(inputs, :yscale,                1.0)
    get!(inputs, :xdisp,                 0.0)
    get!(inputs, :ydisp,                 0.0)
    get!(inputs, :filter_type,           "erf")
    get!(inputs, :lplot_surf3d,          false)
    get!(inputs, :lvolume3d,            false)
    get!(inputs, :smoothing_factor,      0.1)
    get!(inputs, :lrestart,             false)
    get!(inputs, :ndiagnostics_outputs,  0)
    get!(inputs, :Δt,                    0.1)
    get!(inputs, :restart_time,          0.0)
    get!(inputs, :tinit,                 0.0)
    get!(inputs, :tend,                  0.0)
    get!(inputs, :lexact_integration,    false)
    get!(inputs, :llump,                 false)
    get!(inputs, :Nsamp,                 1)
    get!(inputs, :lelementLearning,      false)
    get!(inputs, :ode_adaptive_solver,   false)
    get!(inputs, :output_dir,            "none")
    get!(inputs, :loutput_pert,          false)
    get!(inputs, :lwrite_initial,        false)
    get!(inputs, :μ,                     Float64(0.0))
    get!(inputs, :C1,                    0.0)
    get!(inputs, :C2,                    0.0)
    get!(inputs, :Pr,                    0.7)
    get!(inputs, :lvisc,                 false)
    get!(inputs, :visc_model,            AV())
    get!(inputs, :lrichardson,           false)
    get!(inputs, :lsaturation,           false)
    get!(inputs, :rconst,                Float64(0.0))
    get!(inputs, :iconst,                Int32(1))
    get!(inputs, :lneumann,             false)
    get!(inputs, :penalty,               Float16(0.0))
    get!(inputs, :lsponge,               false)
    get!(inputs, :zsponge,               14000.0)
    get!(inputs, :lmoist,                false)
    get!(inputs, :lprecip,               false)
    get!(inputs, :energy_equation,       "theta")
    get!(inputs, :δtotal_energy,         0.0)
    get!(inputs, :CL,                    CL())
    get!(inputs, :loverwrite_output,     false)
    get!(inputs, :SOL_VARS_TYPE,         TOTAL())
    get!(inputs, :sol_vars_names,        ("q1", "q2", "q3", "q4"))
    get!(inputs, :lsource,               false)
    get!(inputs, :luser_function,        false)
    get!(inputs, :ldss_differentiation,  false)
    get!(inputs, :ldss_laplace,          false)
    get!(inputs, :lamr,                  false)
    get!(inputs, :ladapt,               false)
    get!(inputs, :linitial_refine,       false)
    get!(inputs, :init_refine_lvl,       0)
    get!(inputs, :amr_max_level,         0)
    get!(inputs, :enable_coupling_ready, false)
    get!(inputs, :enable_coupling,       false)
    get!(inputs, :Δt_couple,             :Δt)
    get!(inputs, :couple_time_tol,       1.0e-12)
    get!(inputs, :lwarmup,               false)
end

function mod_inputs_user_inputs!(inputs, rank = 0)

    #Store parsed arguments xxx into inputs[:xxx]
    _parsedToInputs(inputs, parsed_equations, parsed_equations_case_name)

    print_rank(GREEN_FG(string(" # Read inputs dict from ", user_input_file, " ... \n")); msg_rank = rank)
    if rank == 0
        pretty_table(inputs; sortkeys=true, border_crayon = crayon"yellow")
    end

    #
    # Check that necessary inputs exist in the Dict inside .../IO/user_inputs.jl
    #
    mod_inputs_check(inputs, :nop, Int8(4), "w")  #Polynomial order

    # Apply all simple scalar defaults
    _apply_defaults!(inputs)

    # backend-dependent globals (must run after _apply_defaults! sets :backend)
    if inputs[:backend] != CPU()
        global TInt   = Int32
        global TFloat = Float32
        global cpu    = false
    end

    # bulk_fluxes: enabling it forces bdy_fluxes on
    if inputs[:bulk_fluxes] == true && inputs[:bdy_fluxes] == false
        inputs[:bdy_fluxes] = true
    end
    
    
    # lstretch: nested sub-defaults only meaningful when lstretch=true
    if inputs[:lstretch] == true
        get!(inputs, :stretch_type, "powerlaw")
        if inputs[:stretch_type] == "fixed_first"
            get!(inputs, :first_zelement_size, 1.0)
        elseif (inputs[:stretch_type] == "fixed_first_twoblocks_weak"  ||
                inputs[:stretch_type] == "fixed_first_twoblocks_strong" ||
                inputs[:stretch_type] == "fixed_first_twoblocks_strong_weak")
            get!(inputs, :first_zelement_size,       1.0)
            get!(inputs, :max_zelement_size_bottom,  1.0)
            get!(inputs, :zlevel_transition,         1000000000.0)
            get!(inputs, :uniform_zelement_size,     1.0)
            get!(inputs, :max_zelement_size_top,     1.0)
        end
    end
    
    # mount_type / lagnesi: set lagnesi only when :mount_type key is absent
    if !haskey(inputs, :mount_type)
        inputs[:lagnesi] = "agnesi"
    end

    # outformat: string → type conversion (only if user supplied a string value)
    if haskey(inputs, :outformat) && inputs[:outformat] isa AbstractString
        s = lowercase(inputs[:outformat])
        if     s == "png";               inputs[:outformat] = PNG()
        elseif s == "ascii";             inputs[:outformat] = ASCII()
        elseif s == "vtk";               inputs[:outformat] = VTK()
        elseif s == "hdf5" || s == "h5"; inputs[:outformat] = HDF5()
        elseif s == "netcdf";            inputs[:outformat] = NETCDF()
        end
    else
        get!(inputs, :outformat, ASCII())
    end

    # radiation_time_step: depends on :Δt (already defaulted by _apply_defaults!)
    get!(inputs, :radiation_time_step, inputs[:Δt] * 100)

    # diagnostics_at_times: depends on :tend (already defaulted by _apply_defaults!)
    if !haskey(inputs, :diagnostics_at_times)
        inputs[:diagnostics_at_times] = inputs[:tend]
    else
        inputs[:ndiagnostics_outputs] = 0
    end
    # interpolation_nodes: string → type conversion (only if user supplied a string)
    if haskey(inputs, :interpolation_nodes) && inputs[:interpolation_nodes] isa AbstractString
        s = lowercase(inputs[:interpolation_nodes])
        if     s == "llg" || s == "gll" || s == "lgl"; inputs[:interpolation_nodes] = LGL()
        elseif s == "lg"  || s == "gl";                 inputs[:interpolation_nodes] = LG()
        elseif s == "cg"  || s == "gc";                 inputs[:interpolation_nodes] = CG()
        elseif s == "cgl" || s == "gcl";                inputs[:interpolation_nodes] = CGL()
        else
            error("""
                    ERROR in user_inputs.jl --> :interpolation_nodes
                        Chose among: "lgl", "lg", "cg", "cgl"
                  """)
        end
    else
        get!(inputs, :interpolation_nodes, LGL())
    end

    # quadrature_nodes: string → type conversion (only if user supplied a string)
    if haskey(inputs, :quadrature_nodes) && inputs[:quadrature_nodes] isa AbstractString
        s = lowercase(inputs[:quadrature_nodes])
        if     s == "llg" || s == "gll" || s == "lgl"; inputs[:quadrature_nodes] = LGL()
        elseif s == "lg"  || s == "gl";                 inputs[:quadrature_nodes] = LG()
        elseif s == "cg"  || s == "gc";                 inputs[:quadrature_nodes] = CG()
        elseif s == "cgl" || s == "gcl";                inputs[:quadrature_nodes] = CGL()
        else
            error("""
                    ERROR in user_inputs.jl --> :quadrature_nodes
                        Chose among: "lgl", "lg", "cg", "cgl"
                  """)
        end
    else
        get!(inputs, :quadrature_nodes, LGL())
    end

    # ode_solver: warn if not set by user
    # https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/
    if !haskey(inputs, :ode_solver)
        @warn """
                        WARNING in user_inputs.jl --> :ode_solver
                            See usable solvers at
                            https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/
                        SSPRK53 will be used by default.
                            """
        inputs[:ode_solver] = SSPRK54()
    end

    # gmsh_filename_c: derive from gmsh_filename if not set
    if !haskey(inputs, :gmsh_filename_c)
        inputs[:gmsh_filename_c] = get(inputs, :gmsh_filename, "none")
    end

    # Grid entries: required checks depend on lread_gmsh
    if !haskey(inputs, :lread_gmsh) || inputs[:lread_gmsh] == false
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
        mod_inputs_check(inputs, :gmsh_filename,   "e")
        mod_inputs_check(inputs, :gmsh_filename_c, "e")
        mod_inputs_check(inputs, :nsd,  Int8(3), "-")
        mod_inputs_check(inputs, :nelx, Int8(2),  "-")
        mod_inputs_check(inputs, :xmin, Float64(-1.0), "-")
        mod_inputs_check(inputs, :xmax, Float64(+1.0), "-")
        mod_inputs_check(inputs, :nely, Int8(2),  "-")
        mod_inputs_check(inputs, :ymin, Float64(-1.0), "-")
        mod_inputs_check(inputs, :ymax, Float64(+1.0), "-")
        mod_inputs_check(inputs, :nelz, Int8(2),  "-")
        mod_inputs_check(inputs, :zmin, Float64(-1.0), "-")
        mod_inputs_check(inputs, :zmax, Float64(+1.0), "-")
    end #lread_gmsh

    # volume_flux: depends on :lkep
    if inputs[:lkep] == true
        get!(inputs, :volume_flux, "ranocha")
    else
        get!(inputs, :volume_flux, nothing)
    end

    # δvisc: derived from :μ — not user-settable
    inputs[:δvisc] = (inputs[:μ] != 0.0) ? 1.0 : 0.0

    # npx/npy/npz: derived from nelx/nely/nelz
    inputs[:npx] = haskey(inputs, :nelx) ? inputs[:nelx] + 1 : UInt8(2)
    inputs[:npy] = haskey(inputs, :nely) ? inputs[:nely] + 1 : UInt8(2)
    inputs[:npz] = haskey(inputs, :nelz) ? inputs[:nelz] + 1 : UInt8(2)
    if     inputs[:nsd] == 1; inputs[:npy] = UInt8(1); inputs[:npz] = UInt8(1)
    elseif inputs[:nsd] == 2; inputs[:npz] = UInt8(1)
    end

    # :AD validation
    if !haskey(inputs, :AD)
        inputs[:AD] = ContGal()
    elseif inputs[:AD] != ContGal() && inputs[:AD] != FD()
        @mystop(" :AD can only be either ContGal() or FD() at the moment.")
    end

    # :case normalization
    if haskey(inputs, :case)
        inputs[:case] = lowercase(inputs[:case])
    else
        inputs[:case] = ""
    end

    print_rank(GREEN_FG(string(" # Read inputs dict from ", user_input_file, " ... DONE\n")); msg_rank = rank)
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

function mod_io_mkoutdir!(inputs)
    
 #--------------------------------------------------------
    # Create output directory if it doesn't exist:
    #--------------------------------------------------------
    user_defined_output_dir = inputs[:output_dir]

    if inputs[:loverwrite_output]
        outstring = string("output")
    else
        outstring = rank == 0 ? string("output-",  Dates.format(now(), "dduyyyy-HHMMSS")) : ""
        outstring = MPI.bcast(outstring, 0, comm)
    end
    if user_defined_output_dir == "none"
        OUTPUT_DIR = joinpath(case_name_dir, outstring)
        inputs[:output_dir] = OUTPUT_DIR
    else
        OUTPUT_DIR = joinpath(user_defined_output_dir, parsed_equations, parsed_equations_case_name, outstring)
        inputs[:output_dir] = OUTPUT_DIR
    end
    if !isdir(OUTPUT_DIR)
        mkpath(OUTPUT_DIR)
    end

    #--------------------------------------------------------
    # Create restart output/inupt directory if it doesn't exist:
    #--------------------------------------------------------
    if (!haskey(inputs, :restart_output_file_path))
        inputs[:restart_output_file_path] = joinpath(OUTPUT_DIR,string("restart"))
    end

    if (haskey(inputs, :lrestart))
        if(inputs[:lrestart] == true && !haskey(inputs, :restart_input_file_path))
            inputs[:restart_input_file_path] = inputs[:restart_output_file_path]
        end
    else
        inputs[:lrestart] = false
    end

     #--------------------------------------------------------
    # Save a copy of user_inputs.jl for the case being run
    #--------------------------------------------------------
    if rank == 0
        cp(user_input_file, joinpath(OUTPUT_DIR, basename(user_input_file)); force = true)
    end
    
    
    return OUTPUT_DIR 
end
