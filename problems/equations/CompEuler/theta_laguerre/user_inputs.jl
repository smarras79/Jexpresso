function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 1000.0,
        :ode_solver           => SSPRK33(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :tinit                => 0.0,
        :Δt                   => 0.05,
        :ndiagnostics_outputs => 2,
        :case                 => "rtb",
        :lsource              => true, 
        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        #:lexact_integration  => true,
        #:llump               => true,
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        :nop_laguerre        => 16,
        :xfac_laguerre       => 0.0,
        :yfac_laguerre       => 110.0,
        :llaguerre_bc        => true,
        :laguerre_tag        =>"free_slip",
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => [1, 2, 3, 4],
        :μ                   => [0.0, 125.0, 125.0, 125.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB20x20_Lag.msh", #for nop=4
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        :mu_x                => 0.01,
        :mu_y                => 0.01,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :loverwrite_output   => true,
        :outformat           => "vtk",
        :output_dir          => "./output",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
