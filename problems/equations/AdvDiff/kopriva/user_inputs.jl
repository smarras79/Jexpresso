function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 10.0, #2π,
        :Δt                   => 0.005,#8.75e-4,
        :ode_solver           => SSPRK54(),
        :diagnostics_at_times => [0.5, 1, 2, 4],
        :output_dir          => "./output/",
        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :nop                 => 4,      # Polynomial order
        :lsource             => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :μ                    => [0.1, 0.1],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename         => "./meshes/gmsh_grids/kopriva.msh",
        :gmsh_filename         => "./meshes/gmsh_grids/kopriva_periodic.msh",
        #---------------------------------------------------------------------------
        # grid modification parameters
        #--------------------------------------------------------------------------- 
        :xscale              => 1.0,
        :yscale              => 1.0,
        :xdisp               => 0.0,
        :ydisp               => 0.0,
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        #:lwarp               => true,
        #:mount_type          => "agnesi",
        #:a_mount             => 1000.0,
        #:h_mount             => 1.0,
        #:c_mount             => 0.0,
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.15,
        #:mu_y                => 0.15,
        #:filter_type         => "erf",  ##default is erf, use either "erf" for Boyd-Vandeven,"exp" for Warburton Exponential filter, or "quad" for Fischer quadratic filter
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat         => "png",
        :loverwrite_output => true,
        :loutput_pert      => true,  #this is only implemented for VTK for now
        :output_dir        => "./output/",
        #:plot_hlines      => [10.0],
        :loutput_pert      => true
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
