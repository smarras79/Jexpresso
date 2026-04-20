function user_inputs()
    # Solver parameters
    solver_par = Dict(:restart  => true,
                      :memory   => 100,
                      :verbose  => 1,
                      :atol     => 1.e-06,
                      :rtol     => 1.e-06,
                      :itmax    => 1000,
                      )

    # Preconditioner parameters
    prec_sp = Dict(
        :precision    => Float32,
        :ilu_tol      => 0.01,
        :prec_type    => "ilu",
        )

    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => "GMRES", #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.4,
        :tinit                => 0.0,
        :tend                 => 400.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        :restart_input_file_path => "./output/CompEuler/theta/output",
        :ndiagnostics_outputs => 2,
        :case                 => "rtb",
        :lsource              => true, 
        #:backend              => CUDABackend(),
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc               => true, #false by default
        :ivisc_equations     => [1, 2, 3, 4, 5],
        :μ                   => [0.0, 20.0, 20.0, 20.0, 60.0],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x1x1.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x1x10.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x2.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_20x1x20.msh",
        :extra_dimensions    => 1,
        :adaptive_extra_meshes => false,
        :extra_dimensions_order => 4,
        :extra_dimensions_nelemx => 4,
        #:extra_dimensions_nelemy => 4,
        #---------------------------------------------------------------------------
        # Refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => true,
        :init_refine_lvl     => 1,
        #---------------------------------------------------------------------------
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
        :outformat           => "vtk", #"hdf5",
        :output_dir          => "./output/",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
        # Solver parameters
        #---------------------------------------------------------------------------
        :sp                 => solver_par,
        :solver_precision   => Float64,
        :prec_sp            => prec_sp,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
