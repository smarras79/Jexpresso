function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :equation_set => "ns",
        :problem      => "none",
        :tend         => 10,
        :lexact_integration => false,
        :lread_gmsh   => true,
        #:gmsh_filename => "./demo/gmsh_grids/hexa_TFI_2x2.msh",
        :gmsh_filename => "./demo/gmsh_grids/hexa_TFI_10x10.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_TFI_25x25.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_TFI_1x1.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR_coarse.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR_refine.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR_refine_coarse.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock-2x1x1.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock-1x1x1.msh",
        :xmin_bc  => "periodic", #Use either dirichlet or periodic
        :ymin_bc  => "periodic", #Use either dirichlet or periodic
        :zmin_bc  => "periodic", #Use either dirichlet or periodic
        :xmax_bc  => "periodic", #Use either dirichlet or periodic
        :ymax_bc  => "dirichlet", #Use either dirichlet or periodic
        :zmax_bc  => "dirichlet", #Use either dirichlet or periodic
        :bc_exact_xmin => [0.0 0.0 0.0],
        :bc_exact_xmax => [0.0 0.0 0.0],
        :bc_exact_ymin => [0.0 0.0 0.0],
        :bc_exact_ymax => [0.0 0.0 0.0],
        :bc_exact_zmin => [0.0 0.0 0.0],
        :bc_exact_zmax => [0.0 0.0 0.0],
        :nsd                 => 2,           #number of space dimensions
        :interpolation_nodes =>"lgl",        #Choice: lgl, cgl 
        :nop                 => 8,         #Polynomila order
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
