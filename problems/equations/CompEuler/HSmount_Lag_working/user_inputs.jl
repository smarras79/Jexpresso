function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 500.0, #2π,
        :Δt                   => 0.1,#8.75e-4,
        :ode_solver           => SSPRK54(),
        :ndiagnostics_outputs => 100,
        :case                 => "rtb",
        #:CL                   => NCL(),
        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :nop                 => 4,      # Polynomial order
        :nop_laguerre        => 14,
        :yfac_laguerre       => 500.0,
        :luser_bc            => true,
        :lsource             => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => (2,3,4),
        :μ                    => (0.0,1.0,1.0,1.0), #kinematic viscosity constant for θ equation
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/agnesi240kmX30km_coarse.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10_laguerre_top.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/agnesi-120kmx30km-hm5000.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_RTB.msh",
        #:gmsh_filename         => "./meshes/gmsh_grids/hexa_TFI_120x31_top_lag.msh",
        #:gmsh_filename         => "./meshes/gmsh_grids/hexa_TFI_120x31_top_lag_periodic.msh",
        :gmsh_filename         => "./meshes/gmsh_grids/hexa_TFI_120x21_top_lag_periodic.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/agnesi240kmX30km_coarse_laguerreTopLateral.msh",
        #---------------------------------------------------------------------------
        # grid modification parameters
        #---------------------------------------------------------------------------
        :xscale              => 240000.0,
        :yscale              => 15000.0,
        :xdisp               => 0.0,
        :ydisp               => 1.0,
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        :lwarp               => true,
        :mount_type          => "agnesi",
        :a_mount             => 10000.0,
        :h_mount             => 1.0,
        :c_mount             => 0.0,
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => true,
        :mu_x                => 0.2,
        :mu_y                => 0.2,
        :filter_type         => "erf",  ##default is erf, use either "erf" for Boyd-Vandeven,"exp" for Warburton Exponential filter, or "quad" for Fischer quadratic filter
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
