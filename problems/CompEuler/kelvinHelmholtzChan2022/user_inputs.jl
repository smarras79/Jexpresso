function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
	:ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 1.0e-2,
        :tinit                => 0.0,
        :tend                 => 3.0,
        :diagnostics_at_times => (0.0:1.0:3.0),
        :restart_time         => 0.0,
        :lrestart             => false,
        #:case                 => "rtb",
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(), #THETA(),
        :use_named_tuples     => true, # Converts inputs to named tuples
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 3,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :μ                   => [0.0, 1.0, 1.0, 1.0], #horizontal viscosity constant for momentum
        #:μ                   => [0.0, 0.25, 0.25, 0.25], #horizontal viscosity constant for momentum
        :visc_model           => SMAG(),
        #:visc_model           => VREM(),
        #---------------------------------------------------------------------------
        # LKEP:
        #---------------------------------------------------------------------------
        :lkep        => true,
	:entropy_variables => false,
        :volume_flux => kennedy_gruber(),
        #:volume_flux => artiano_tec(),
        #:volume_flux => central_theta(),
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
        :lwrite_initial      => false,
       # :lwrite_initial      => true,
        :output_dir          => "./output-theta/",
        #:output_dir          => "./test/CI-run",
        :loutput_pert        => false,  #this is only implemented for VTK for now
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
