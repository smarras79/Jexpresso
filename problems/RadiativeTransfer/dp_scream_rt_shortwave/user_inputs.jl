function user_inputs()
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
        :gmsh_filename       => "./JexpressoMeshes/meshes/gmsh_grids/RICO_5x5.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_20x1x20.msh",
        :lRT_problem        => true,
        :lRT_from_data       => true,
        :RT_data_file       => "./DP_SCREAM_DATA/RICO_5x5/scream_dpxx_RICO_5x5/test/input/scream_dpxx_RICO_5x5.scream.INSTANT.nhours_x1.2004-12-16-00000.012.in.nc",
        :RT_shortwave   => true,
        :RT_longwave    => false,
        :RT_μ0 => cos(deg2rad(83.0)),
        :RT_ϕ0 => 0.0,
        :RT_radiative_heating => true,
        :extra_dimensions    => 2,
        :adaptive_extra_meshes => false,
        :RT_amr_threshold => [0.99999],
        :extra_dimensions_order => 3,
        :extra_dimensions_nelemx => 3,
        :extra_dimensions_nelemy => 3,
        :rad_HG_g => 0.85,
        :lcubed_sphere_angular_mesh => false,
        :extra_dimensions_xmax => π,
        :extra_dimensions_ymax => 2*π,
        :extra_dimensions_xmin => 0,
        #:extra_dimensions_nelemy => 4,
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
        :outformat           => "netcdf", #"hdf5",
        :output_dir          => "./output/",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
