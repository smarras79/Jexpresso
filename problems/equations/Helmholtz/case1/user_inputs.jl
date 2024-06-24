function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 1000.0,
        :ode_solver           => "BICGSTABLE", #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :ndiagnostics_outputs => 1,
        :lsource              => true, 
        :llinsolve            => true,
        #:backend              => MetalBackend(),
        #:CL                   => NCL(), #CL() is defaults
        #:SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        #:lexact_integration  => true,
        #:llump               => true,
        :interpolation_nodes =>"lgl",
        :nop                 => 10,      # Polynomial order
        :nop_laguerre        => 14,
        :xfac_laguerre       => 0.25,
        :yfac_laguerre       => 0.0,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lvisc                => true, #false by default NOTICE: works only for Inexact
        #:ivisc_equations      => (1, 2, 3, 4),
        #:μ                   => (0.0, 75.0, 75.0, 75.0), #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB_unitsize.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_helmholtz.msh", #for nop=4
        #:gmsh_filename       => "./meshes/gmsh_grids/mygmsh.msh", #for nop=4
        #:gmsh_filename        => "./meshes/gmsh_grids/test_allocation.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_120x31_periodic.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB27x27.msh", #for nop=3
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB40x40.msh", #for nop=2
        #---------------------------------------------------------------------------
        # grid modification parameters
        #---------------------------------------------------------------------------
        :xscale              => 5.0,
        :yscale              => 3.14,
        :xdisp               => 1.0,
        :ydisp               => 0.0,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "png",
        :output_dir          => "./output/",
        :plot_vlines         => [5.0],
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
