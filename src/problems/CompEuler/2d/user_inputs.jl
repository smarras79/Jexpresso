function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 2.0, #2π,
        :Δt                   => 1.0e-4,#8.75e-4,
        :ode_solver           => "SSPRK53",
        :ndiagnostics_outputs => 100,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :lexact_integration  => false,
        :nop                 => 5,      # Polynomial order
        :luser_bc            => true,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        #:lplot_surf3d        => true,   #false by default
        #:smoothing_factor    => 1.0, #factor for spline2d interpolation. 
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #true by default
        :νx                   => 0.00, #kinematic viscosity constant
        :νy                   => 0.00, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_1x1.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x2.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_44x20.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_32x16.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_25x25.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_UNSTR_coarse.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10_periodic.msh",
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
