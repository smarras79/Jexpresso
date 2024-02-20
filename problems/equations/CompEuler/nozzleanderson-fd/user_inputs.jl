function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 0.1,
        :ode_solver           => ABM54(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :Δt                   => 10e-6,
        :ndiagnostics_outputs => 20,
        :case                 => "rtb",
        :lsource              => true, 
        #:CL                   => NCL(), #CL() is defaults
        #:SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        #:lexact_integration  => true,
        #:llump               => true,
        :interpolation_nodes =>"lgl",
        :nop                 => 1,
        :AD                  => FD(),
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lvisc                => true, #false by default NOTICE: works only for Inexact
        #:νx                   => 30.0, #horizontal viscosity constant for momentum
        #:νy                   => 30.0, #vertical   viscosity constant for momentum
        #:κ                    => 60.0, #kinematic viscosity constant for θ equation
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        #:lread_gmsh => false, #If false, a 1D problem will be enforced
        :xmin => 0.0,
        :xmax => 3.0,
        :nelx => 100,
        :lperiodic_1d => false,
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
        :outformat  => "png",
        :output_dir => "./output/",
        :outvars    => ("ρ", "u", "T", "p"),
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end