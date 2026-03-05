function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 1.0,
        :ode_solver           => "BICGSTABLE", #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :ndiagnostics_outputs => 1,
        :lsource              => true, 
        :llinsolve            => true,
        :ldss_laplace         => true,
        :lelementLearning     => true,
        :lEL_Train            => true,
        :Nsamp                => 1,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_1x1.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_3x3.msh",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 12,      # Polynomial order
        #:nop_laguerre        => 14,
        #:xfac_laguerre       => 0.25,
        #:yfac_laguerre       => 0.0,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lvisc                => true, #false by default NOTICE: works only for Inexact
        #:ivisc_equations      => (1, 2, 3, 4),
        #:μ                   => (0.0, 75.0, 75.0, 75.0), #horizontal viscosity constant for momentum
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
        :outformat           => "vtk",
        :output_dir          => "./output/",
        :loverwrite_output   => true,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
