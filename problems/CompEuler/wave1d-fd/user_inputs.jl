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
        :ndiagnostics_outputs => 10, #these are steps, not seconds
        :output_dir          => "./",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :AD                  => FD(),
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 1,     # Polynomial order
        :lexact_integration  => false,
        :lsource             => false,
        :lperiodic_1d        => true, #false by default
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => 1,
        :μ                   => 0.1, #horizontal viscosity constant for momentumx
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
        :xmax          =>   5.0,
        :nelx          =>   50,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
