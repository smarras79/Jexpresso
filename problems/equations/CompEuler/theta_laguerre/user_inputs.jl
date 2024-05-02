function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 500.0,
        :ode_solver           => SSPRK33(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :tinit                => 0.0,
        :Δt                   => 0.05,
        :ndiagnostics_outputs => 20,
        :case                 => "rtb",
        # :ndiagnostics_outputs => 6,
        :lsource              => true, 
        #:CL                   => NCL(), #CL() is defaults
        :backend              => MetalBackend(),
        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        #:lexact_integration  => true,
        #:llump               => true,
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        :nop_laguerre        => 14,
        :xfac_laguerre       => 0.0,
        :yfac_laguerre       => 110.0,
        :llaguerre_bc        => true,
        :laguerre_tag        =>"free_slip",
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => (1, 2, 3, 4),
        :μ                   => [0.0, 75.0, 75.0, 75.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB_unitsize.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB20x20_Lag.msh", #for nop=4
        #:gmsh_filename       => "./meshes/gmsh_grids/mygmsh.msh", #for nop=4
        #:gmsh_filename        => "./meshes/gmsh_grids/test_allocation.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_120x31_periodic.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB27x27.msh", #for nop=3
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB40x40.msh", #for nop=2
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
        :output_dir          => "./output/",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
