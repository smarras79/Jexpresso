function user_inputs()
    
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.0005,
        :tinit                => 0.0,
        :tend                 => 20.0,
        :diagnostics_at_times => range(0,20,41),
        :case                 => "shu",
        #:case                 => "hw",
        #:case                 => "vincent",
        :lsource              => true,
        :lperiodic            => true,
        #:backend              => MetalBackend(),
        #:SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => [1, 2, 3, 4],
        :μ                    => [0.0, 0.05, 0.05, 0.05], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/vortex_20x20.msh", #for nop=4
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.01,
        #:mu_y                => 0.01,
        #:filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "png",
        :loverwrite_output   => true,
        :lwrite_initial      => true,
        :output_dir          => "./output",
        #:output_dir          => "./test/CI-run",
        #:loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
