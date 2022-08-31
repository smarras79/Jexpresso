function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :equation_set => "ns",
        :problem      => "rtb",
        :lexact_integration => false,
        :lread_gmsh   => false,
        #:gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR_coarse.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock-2x1x1.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock-1x1x1.msh",
        :nsd          => 1,           #number of space dimensions
        :nop          => 1,           #Polynomila order
        :nelx         => 1,           #N. elements in x
        :nely         => 0,           #N. elements in y
        :nelz         => 0,           #N. elements in z
        :xmin         => 0,
        :xmax         => 5,
        :ymin         =>-1,
        :ymax         => 1,
        :zmin         =>-1,
        :zmax         => 1
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
