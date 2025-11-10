function user_inputs()
    
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
<<<<<<< HEAD
	:ode_solver           => SSPRK43(), #RDPK3SpFSAL49(),#CarpenterKennedy2N54(),#SSPRK54(),  #ORK256(), #SSPRK54(),
	:ode_adaptive_solver  => true,
        :Δt                   => 0.005,
        :tinit                => 0.0,
        :tend                 => 20.0,
        :diagnostics_at_times => (0.0:0.5:20.0),
=======
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.0003,
        :tinit                => 0.0,
        :tend                 => 10.0,
        :diagnostics_at_times => (0.0:0.25:10.0),
>>>>>>> bd48c1fbb64854411f292b3294ac887d132e8dc6
        :restart_time         => 0.0,
        :lrestart             => false,
        :restart_input_file_path => "/home/leon/njit/Jexpresso_gigales/Jexpresso/problems/equations/CompEuler/theta",
        :case                 => "rtb",
<<<<<<< HEAD
        :lsource              => false,
	#:SOL_VARS_TYPE        => THETA(), #PERT(), #TOTAL() is default
=======
        :lsource              => false, 
        :SOL_VARS_TYPE        => THETA(),
>>>>>>> bd48c1fbb64854411f292b3294ac887d132e8dc6
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
<<<<<<< HEAD
        :nop                 => 6,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lkep                => true,
        :lvisc               => true, #false by default NOTICE: works only for Inexact
=======
        :nop                 => 7,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
>>>>>>> bd48c1fbb64854411f292b3294ac887d132e8dc6
        :μ                   => [0.0, 1.0, 1.0, 1.0], #horizontal viscosity constant for momentum
        #:μ                   => [0.0, 0.25, 0.25, 0.25], #horizontal viscosity constant for momentum
        :visc_model           => SMAG(),
        #---------------------------------------------------------------------------
        # LKEP:
        #---------------------------------------------------------------------------
        #:lkep        => true,
        :volume_flux => "ec",
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10_unitsquare.msh", #for nop=4
	#:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_40x40_unitsquare.msh", #for nop=4
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_32x32_unitsquare.msh", #for nop=4
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
        :loverwrite_output   => false,
        :lwrite_initial      => true,
<<<<<<< HEAD
        :output_dir          => "./parametric_runs/nop6_flux_kennedy_gruber_test",
        #:output_dir          => "./parametric_runs/nop7_flux_kg",
=======
        :output_dir          => "./output-theta/",
        #:output_dir          => "./test/CI-run",
>>>>>>> bd48c1fbb64854411f292b3294ac887d132e8dc6
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
