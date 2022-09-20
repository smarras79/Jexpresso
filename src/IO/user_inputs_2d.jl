function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :equation_set => "ns",
        :problem      => "none",
        :tend         => 2.5,
        :lexact_integration => false,
        #:lread_gmsh   => false,
        #:gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR.msh",
        :gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR_coarse.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock-2x1x1.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock-1x1x1.msh",
        :nsd                 => 1,           #number of space dimensions
        :interpolation_nodes =>"lgl",        #Choice: lgl, cgl 
        :nop                 => 100,         #Polynomila order
        :nelx                => 1,           #N. elements in x
        :nely                => 1,           #N. elements in y
        :nelz                => 1,           #N. elements in z
        :xmin                => -1,
        :xmax                => 1,
        :ymin                => 0,
        :ymax                => 0,
        :zmin                => 0,
        :zmax                => 0
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
