function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :equation_set => "burgers",
        :problem      => "burgers1d",
        :nsd          => 1,  #number of space dimensions
        :npx          => 100,  #N. points in x
        :npy          => 1,  #N. points in y
        :npz          => 1,  #N. points in z
        :xmin         =>0,
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
