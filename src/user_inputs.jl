function user_inputs()
    inputs = (
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        equation_set = "burgers",
        problem      = "burgers1d",
        nsd          = 1,   #number of space dimensions
        npx          = 100, #N. points in x
        npy          = 100, #N. points in y
        npz          = 100, #N. points in z
        xmin         =  -1,
        xmax         =   1,
        ymin         =  -1,
        ymax         =   1,
        zmin         =  -1,
        zmax         =   1
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------
    
    return inputs
    
end
