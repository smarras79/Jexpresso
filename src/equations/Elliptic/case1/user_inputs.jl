function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :ode_solver          => "BICGSTABLE",
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 4,     # Polynomial order
        :lexact_integration  => false,
        :luser_bc            => false,
        #---------------------------------------------------------------------------
        # Output parameters:
        #---------------------------------------------------------------------------
        :outformat           => "png",
        #:lplot_surf3d        => true, #default is false. when true, set a value for :smoothing_factor. as small as possible but not too small!
        :smoothing_factor    => 1.0e-4, #factor for spline2d interpolation. 
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:νx                   => 0.01, #kinematic viscosity constant
        #:νy                   => 0.01, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/plate_hole.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/2d-grid-zeroall.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/circle_TFI.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/circle1.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x2.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_25x25.msh",
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => faluse): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin          =>   0.0,
        :xmax          =>   1.0,
        :nelx          =>   25,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
