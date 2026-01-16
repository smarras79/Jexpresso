function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # User define your inputs below: the order doesn't matter
        # IMPORTANT NOTICE: DO NOT FORGET the "," at the end of each entry!!!
        #---------------------------------------------------------------------------
        :ode_solver          => SSPRK54(),
        :tend                 => 10.0,
        :Δt                   => 0.002,
        :output_dir          => "./",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :interpolation_nodes_extra => "lgl",

        :nop                 => 6,     # Polynomial order
        :nop_extra           => 6,

        :lexact_integration  => false,
        :lexact_integration_extra => false,

        :lsource             => false,

        :lperiodic_1d        => true, 
        :lperiodic_1d_extra  => false,
        
        :luser_bc             => true,
        :l_extra_coord       => true,

        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => false,
        :nsd                 => 1,
        :nelx                => 10,
        :xmin                => 0,
        :xmax                => 10*pi,

        :nsd_extra                 => 1,
        :nelx_extra                => 10,
        :xmin_extra                => -6,
        :xmax_extra                => 6

    ) #Dict

    return inputs
    
end
