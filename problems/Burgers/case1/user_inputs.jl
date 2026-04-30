function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 1D viscous Burgers equation (conservation form):
        #
        #     ∂q/∂t + ∂F(q)/∂x = ν ∂²q/∂x²,    F(q) = q²/2
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK33(),
        :tend                 => 1.0,
        :Δt                   => 1.0e-4,
        :diagnostics_at_times => (0:0.1:1.0),
        :output_dir          => "./",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 4,     # Polynomial order
        :lexact_integration  => false,
        :lsource             => false,
        :lperiodic_1d        => true, #false by default
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc               => true,     #turn on physical viscosity
        :μ                   => [0.01],   #kinematic viscosity ν for the Burgers equation
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => false, #If false, a 1D problem will be enforced
        #---------------------------------------------------------------------------
        # Output formats: "png" -> plots to png file. "ascii" -> data to npoin file
        #---------------------------------------------------------------------------
        :outformat           => "png",
        :loverwrite_output   => true,
        :output_dir          => "./output",
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => false): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin          =>   0.0,
        :xmax          =>   1.0,
        :nelx          =>   50,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs

end
