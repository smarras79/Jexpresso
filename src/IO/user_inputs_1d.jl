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
        :nsd                 => 1,    #number of space dimensions
        :interpolation_nodes =>"lgl", #Choice: lgl, cgl 
        :nop                 => 8,    #Polynomila order
        #----------------------------------------------
        # Build native 1D grid.
        # For 2D/3D read a GMSH grid instead
        #----------------------------------------------
        :nelx                => 1,    #N. elements in x
        :xmin                => -1,
        :xmax                => 1,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
