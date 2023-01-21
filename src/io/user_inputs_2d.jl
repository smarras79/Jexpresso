function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :equation_set        => "NS",
        :problem             => "none",
        :tend                => π,
        :Δt                  => 0.001,
        :nsd                 => 2,      #number of space dims is overwritten by an external grid
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :lexact_integration  => false,
        :nop                 => 5,      # Polynomila order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :νx                   => 0.1, #kinematic viscosity constant
        :νy                   => 0.1, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_TFI_2x2.msh",
        :gmsh_filename       => "./demo/gmsh_grids/hexa_TFI_10x10.msh",
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_TFI_25x25.msh",
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_TFI_1x1.msh",
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_oneblock.msh",
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_UNSTR_coarse.msh",
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_UNSTR_refine.msh",
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_UNSTR_refine_coarse.msh",
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_oneblock-2x1x1.msh",
        #:gmsh_filename       => "./demo/gmsh_grids/hexa_oneblock-1x1x1.msh",
        #---------------------------------------------------------------------------
        # Boundary conditions:
        #---------------------------------------------------------------------------
        :xmin_bc       => "periodic", #Use either dirichlet or periodic
        :ymin_bc       => "periodic", #Use either dirichlet or periodic
        :zmin_bc       => "periodic", #Use either dirichlet or periodic
        :xmax_bc       => "periodic", #Use either dirichlet or periodic
        :ymax_bc       => "periodic", #Use either dirichlet or periodic
        :zmax_bc       => "periodic", #Use either dirichlet or periodic
        :bc_exact_xmin => [0.0 0.0 0.0],
        :bc_exact_xmax => [0.0 0.0 0.0],
        :bc_exact_ymin => [0.0 0.0 0.0],
        :bc_exact_ymax => [0.0 0.0 0.0],
        :bc_exact_zmin => [0.0 0.0 0.0],
        :bc_exact_zmax => [0.0 0.0 0.0],
        
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
