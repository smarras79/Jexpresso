function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        :ode_solver          => SSPRK53(),
        :tend                 => 3.0,
        :Δt                   => 1.0e-4,
        #:ndiagnostics_outputs => 30, #these are steps, not seconds
        :diagnostics_at_times => range(0.0, 3.0, 60),
        :output_dir          => "./",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 4,     # Polynomial order
        :lsource             => false,
        :lperiodic_1d        => false, #false by default
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
#        :lvisc               => true,
        :ivisc_equations     => [1, 2, 3],
        :μ                   => [0.0, 0.01, 0.01],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => false, #If false, a 1D problem will be enforced
        #---------------------------------------------------------------------------
        # Output formats: "png" -> plots to png file. "ascii" -> data to npoin file
        #---------------------------------------------------------------------------
        :outformat         => "png", #choice: "png", "ascii" (default is ascii)
        :loverwrite_output => true,
        #:output_dir        => "./test/CI-ref/", 
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => faluse): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin          => -5.0,
        :xmax          =>  5.0,
        :nelx          => 50,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
