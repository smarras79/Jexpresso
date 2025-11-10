function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
	:ode_solver           => SSPRK43(), #RDPK3SpFSAL49(),#CarpenterKennedy2N54(),#SSPRK54(),  #ORK256(), #SSPRK54(),
	:ode_adaptive_solver  => true,
        :Δt                   => 0.005,
        :tinit                => 0.0,
        :tend                 => 10.0,
        :diagnostics_at_times => (0.0:0.5:10.0),
        :restart_time         => 0.0,
        :lrestart             => false,
        :restart_input_file_path => "/home/leon/njit/Jexpresso_gigales/Jexpresso/problems/equations/CompEuler/theta",
        :case                 => "rtb",
        :lsource              => false,
	:SOL_VARS_TYPE        => THETA(), #PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
       :lkep                => true,
        :lvisc               => true, #false by default NOTICE: works only for Inexact
        :μ                   => [0.0, 0.25, 0.25, 0.25], #horizontal viscosity constant for momentum
        #:μ                   => [0.0, 0.0001, 0.0001, 0.0001], #horizontal viscosity constant for momentum
        #:visc_model           => VREM(),
        :visc_model    => SMAG(), #AV(), #SMAG(), AV(), DSMAG(), VREM()
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_20x20_unitsquare.msh", #for nop=4
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10_unitsquare.msh", #for nop=4
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_40x40_unitsquare.msh", #for nop=4
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_32x32_unitsquare.msh", #for nop=4
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => false,
        :mu_x                => 0.05,
        :mu_y                => 0.05,
        :filter_type         => "erf",
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
        :init_refine_lvl     => 2,
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
