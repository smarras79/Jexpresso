function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Stage 0 Test Case: RT with Spatial AMR Detection
        # This configuration tests the spatial non-conformity detection added in Stage 0
        #---------------------------------------------------------------------------
        :ode_solver           => "GMRES",
        :Δt                   => 0.4,
        :tinit                => 0.0,
        :tend                 => 0.4,  # Very short - just for detection test
        :ndiagnostics_outputs => 1,
        :case                 => "rtb",
        :lsource              => true,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Low polynomial order for fast testing
        #---------------------------------------------------------------------------
        # Physical parameters
        #---------------------------------------------------------------------------
        :lvisc               => true,
        :ivisc_equations     => [1, 2, 3, 4, 5],
        :μ                   => [0.0, 20.0, 20.0, 20.0, 60.0],
        #---------------------------------------------------------------------------
        # Mesh parameters: Use small uniform mesh for testing
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./JexpressoMeshes/meshes/gmsh_grids/hexa_TFI_3d_rad.msh",
        #---------------------------------------------------------------------------
        # Extra dimensions (angular mesh)
        #---------------------------------------------------------------------------
        :extra_dimensions    => 2,
        :lRT_problem        => true,
        :lmanufactured_solution => true,
        :adaptive_extra_meshes => false,  # No angular adaptation for this test
        :extra_dimensions_order => 2,
        :extra_dimensions_nelemx => 2,
        :extra_dimensions_nelemy => 2,
        :lcubed_sphere_angular_mesh => false,
        :extra_dimensions_xmax => π,
        :extra_dimensions_ymax => 2*π,
        :extra_dimensions_xmin => 0,
        #---------------------------------------------------------------------------
        # STAGE 0 TEST CONFIGURATION: Spatial AMR Detection
        #---------------------------------------------------------------------------
        # Enable spatial refinement in driver for testing
        :lRT_spatial_amr     => true,   # Driver will selectively refine mesh
        :ladapt              => true,   # Enable AMR framework from start
        :linitial_refine     => false,  # Not using uniform initial refinement
        :lamr                => true, 
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output/RT_Stage0_Test/",
        :loutput_pert        => true,
    )

    return inputs
end
