function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 2.0, #2π,
        :Δt                   => 5.0e-4,#8.75e-4,
        :ode_solver           => "SSPRK53",
        :ndiagnostics_outputs => 2,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :lexact_integration  => false,
        :nop                 => 4,      # Polynomial order
        :luser_bc            => true,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "png",
        #:lplot_surf3d        => true,   #false by default
        #:smoothing_factor    => 1.0, #factor for spline2d interpolation. 
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, # default is true
        :νx                   => 0.0, #kinematic viscosity constant
        :νy                   => 0.0, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :nelx                 =>  50,
        :xmin                 =>  -5.0,
        :xmax                 =>   5.0,
        #:lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_1x1.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x2.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_25x25.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_UNSTR_coarse.msh",
        :gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10_periodic.msh",
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
