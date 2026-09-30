function user_inputs()
    
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(),
        :Δt                   => 0.5,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        :diagnostics_at_times => (0:100:1000),
        :restart_time         => 50000,
        :lrestart             => false,
        #:CL                   => NCL(),
        :restart_input_file_path => "/home/leon/njit/Jexpresso_gigales/Jexpresso/problems/equations/CompEuler/theta",
        :case                 => "rtb",
        :lsource              => true, 
        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc          => true, #false by default NOTICE: works only for Inexact       
        :visc_model     => AV(),
        #:visc_model     => VREM(),
        #:visc_model     => SMAG(),
        :energy_equation => "theta",
        #:μ              => [0.0, 1.0, 1.0, 2.0], #horizontal viscosity constant for momentum
        :μ              => [0.0, 125.0, 125.0, 125.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        #---------------------------------------------------------------------------
        # Built-in Cartesian grid: no GMSH file. The same 10×10 quads on
        # [-5000,5000]×[0,10000] as test/CI-runs/CompEuler/theta/hexa_TFI_10x10.msh,
        # with the same boundary names (the physical names of that file); the
        # solution agrees with the GMSH run to round-off (~1e-11 relative
        # after 2000 steps). See jx_cartesian_model in src/kernel/mesh/mesh.jl.
        #---------------------------------------------------------------------------
        :lcartesian_grid     => true,
        :nelx                => 10, :xmin => -5000.0, :xmax =>  5000.0,
        :nely                => 10, :ymin =>     0.0, :ymax => 10000.0,
        :cartesian_bdy       => Dict(:xmin => "free_slipz", :xmax => "free_slipz",
                                     :ymin => "free_slipx", :ymax => "free_slipx"),
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.01,
        #:mu_y                => 0.01,
        #:filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        :lwrite_initial      => true,
        :output_dir          => "./output",
        #:output_dir          => "./test/CI-run",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :ladapt              => false,
        #---------------------------------------------------------------------------
        # AMR parameters
        #---------------------------------------------------------------------------
        :amr_freq            => 200,
        :amr_max_level       => 2,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
