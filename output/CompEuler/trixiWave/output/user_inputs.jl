function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        :ode_solver          => SSPRK54(),
        :tend                 => 0.4,
        :Δt                   => 1.6781e-04,
        :diagnostics_at_times => (0:0.05:0.4),
        :output_dir          => "./",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 3,     # Polynomial order
        :lsource             => false,
        :lperiodic_1d        => true, #false by default
        #---------------------------------------------------------------------------
        # Entropy/energy preserving discretizations:
        #---------------------------------------------------------------------------
        :lkep => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lvisc               => true,
        :μ                   => [0.001, 0.001, 0.001],
        #:μ                   => [0.0001, 0.0001, 0.0001],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => false, #If false, a 1D problem will be enforced
        #---------------------------------------------------------------------------
        # Output formats: "png" -> plots to png file. "ascii" -> data to npoin file
        #---------------------------------------------------------------------------
        :outformat         => "png", #choice: "png", "ascii" (default is ascii)
        :loverwrite_output => true,
        :output_dir        => "./output",
        #:output_dir        => "./test/CI-ref/",
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => faluse): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin          => 0.0,
        :xmax          => 2.0,
        :nelx          =>  12,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs

end
