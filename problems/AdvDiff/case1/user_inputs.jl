function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        :ode_solver          => SSPRK33(),
        :tend                 => 3.0,
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => (0:0.1:3),
        :output_dir          => "./",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 4,     # Polynomial order
        :lexact_integration  => false,
        :lsource             => false,
        :lperiodic_1d        => true, #false by default
        #:backend             => MetalBackend(),
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false,
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
        :outformat           => "png",
        :loverwrite_output   => true,
        :output_dir          => "./output",
        #:output_dir          => "./test/CI-run",
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => faluse): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin          =>   0.0,
        :xmax          =>   2.0,
        :nelx          =>   50,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
