function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        :ode_solver          => SSPRK33(),
        :tend                 => 9.0,
        :Δt                   => 1.0e-3,
        :ndiagnostics_outputs => 30, #these are steps, not seconds
        :output_dir          => "./output/",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 6,     # Polynomial order
        :nop_laguerre        => 50,
        :lexact_integration  => false,
        :lsource             => true,
        :llaguerre_1d_right  => true,
        :llaguerre_1d_left   => true,
        :laguerre_beta       => 1.0,
        :yfac_laguerre       => 0.05,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false,
        :νx                   => 0.01, #kinematic viscosity constant
        :νy                   => 0.01, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => false, #If false, a 1D problem will be enforced
        #---------------------------------------------------------------------------
        # Output formats: "png" -> plots to png file. "ascii" -> data to npoin file
        #---------------------------------------------------------------------------
        :outformat         => "hdf5", #choice: "png", "ascii" (default is ascii)
        :loverwrite_output => true,
        :output_dir        => "./CI-runs", #this is in github 
        :plot_vlines       => [-2.5,2.5],
        :plot_axis         => [-0.05,0.55, -0.35,0.35],
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => faluse): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin          =>   -2.5,
        :xmax          =>   2.5,
        :nelx          =>   50,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
