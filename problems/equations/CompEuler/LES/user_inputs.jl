function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.04,
        :tinit                => 0.0,
        :tend                 => 10000,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        #:restart_input_file_path => "",
        :diagnostics_at_times => (0:1.0:1000),
        :lsource              => true,
        :sounding_file        => "./data_files/input_sounding_teamx_u10_flat_noheader.dat",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 3,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lwall_model          => true,
        :bdy_fluxes           => true,
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        #:visc_model           => SMAG(),
        :visc_model           => AV(),
        :ivisc_equations      => [1, 2, 3, 4, 5],
        # smagorinsky, cs = 0.23, input cs^2 for momentum cs^2/Pr for other equations, where Pr = 1/3
        #:μ                    => [0.0, 0.53, 0.53, 0.53, 1.6], #horizontal viscosity constant for momentum
        :μ                    => [0.0, 1, 1, 1, 1.6], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/LESICP_10240x10240x5000.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/LESICP_16x1x16.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/LESICP_32x32x16.msh",
        #:linitial_refine     => false,
        :init_refine_lvl     => 1,
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
        :output_dir          => "./output/",
        :loverwrite_output   => true,  #this is only implemented for VTK for now
        :lwrite_initial      => true,
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
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
