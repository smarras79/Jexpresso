function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # CI version: elliptic direct solve — tend is not used for time-stepping
        # Uses smaller mesh for faster CI execution
        #---------------------------------------------------------------------------
        :tend                 => 1.0,
        :ode_solver           => "BICGSTABLE",
        :ndiagnostics_outputs => 1,
        :lsource              => true,
        :llinsolve            => true,
        :lsparse              => true,
        :ldss_laplace         => true,
        :ldss_differentiation => false,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_15x15.msh",
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
        :outformat           => "hdf5",
        :output_dir          => "none",
        :loverwrite_output   => true,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
