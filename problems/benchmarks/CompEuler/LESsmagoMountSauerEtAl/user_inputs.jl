function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.25, #0.065,
        :tinit                => 0.0,
        :tend                 => 7200,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        #:restart_input_file_path => "",
        :diagnostics_at_times => (0:5.0:7200),
        :lsource              => true,
        :sounding_file        => "./data_files/input_sounding_teamx_u10_flat_noheader.dat",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lwall_model          => true,
        :bdy_fluxes           => true,
        :lvisc                => true, #false by default
        :visc_model           => SMAG(),
        #:visc_model           => AV(),
        :ivisc_equations      => [1, 2, 3, 4, 5],
        # smagorinsky, cs = 0.23, input cs^2 for momentum cs^2/Pr for other equations, where Pr = 1/3
        #:μ                    => [0.0, 0.53, 0.53, 0.53, 1.6], #horizontal viscosity constant for momentum
        :μ                    => [0.0, 10, 10, 10, 15], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        :lwarp               => true,
        :mount_type          => "sauer", #"agnesi",
        :a_mount             => 36000.0,
        :h_mount             =>  900.0,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lwarmup          => false,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
        #:gmsh_filename_c  => "./meshes/gmsh_grids/LESICP_32x2x24_zmax3000.msh",
        :gmsh_filename    => "./meshes/gmsh_grids/LESICPsauer.msh",
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
        #:output_dir          => "/scratch/smarras/smarras/output/64x64x36_5kmX5kmX3km/",
	#:output_dir          => "/scratch/smarras/smarras/output/64x64x36_5kmX5kmX3km_MORECORES/",
        #:output_dir          => "/scratch/smarras/smarras/output/64x64x24fewcores/",
        :output_dir          => "./outputSauer/",
        :loverwrite_output   => true,  #this is only implemented for VTK for now
        :lwrite_initial      => true,
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        #:linitial_refine     => true,
        :init_refine_lvl     => 1,
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :ladapt              => false,
        :amr                 => true,
        #---------------------------------------------------------------------------
        # AMR parameters
        #---------------------------------------------------------------------------
        :amr_freq            => 20,
        :amr_max_level       => 1,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
