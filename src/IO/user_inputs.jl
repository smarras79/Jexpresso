function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :equation_set => "ns",
        :problem      => "rtb",
        :lread_gmsh   => false,
        :gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_UNSTR_coarse.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock-2x1x1.msh",
        #:gmsh_filename => "./demo/gmsh_grids/hexa_oneblock-1x1x1.msh",
        :nsd          => 3,           #number of space dimensions
        :nop          => 3,           #Polynomila order        
        :npx          => 10,         #N. points in x
        :npy          => 10,           #N. points in y
        :npz          => 5,           #N. points in z
        :xmin         => 0,
        :xmax         => 2π,
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
