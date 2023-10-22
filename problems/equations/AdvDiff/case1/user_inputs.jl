function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        :tend                 => 2.5,
        :Δt                   => 1.0e-3,
        :ndiagnostics_outputs => 15, #these are steps, not seconds
        :lsource              => false,
        #:output_dir          => "/Users/simone/runs/",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :ode_solver          => "AB4", #"Tsit5",
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 4,     # Polynomial order
        :lexact_integration  => false,
        :luser_bc            => true,
        :lsource             => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lvisc                => false,
        :νx                   => 0.01, #kinematic viscosity constant
        :νy                   => 0.01, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => false, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/2d-grid.msh", 
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_25x25.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/circle_TFI.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10_periodic.msh",
        #---------------------------------------------------------------------------
        # Output formats: "png" -> plots to png file. "ascii" -> data to npoin file
        #---------------------------------------------------------------------------
        :outformat     => "png", #choice: "png", "ascii" (default is ascii)
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => faluse): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin          =>   0.0,
        :xmax          =>   2.0,
        :nelx          =>   50,
        #:lperiodic1d   => true,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
