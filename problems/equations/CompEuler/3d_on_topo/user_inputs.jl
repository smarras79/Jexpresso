function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.4,
        :tinit                => 0.0,
        :tend                 => 0.8,
        :diagnostics_at_times => [0.4, 0.8],
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        :restart_input_file_path => "./output",
        :case                 => "rtb",
        :lsource              => true, 
        #:backend              => MetalBackend(),
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => [1, 2, 3, 4, 5],
        :μ                   => [0.0, 80.0, 80.0, 80.0, 240.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x1x1.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x1x10.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10x10.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_20x1x20.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_DOMEX_coarse.msh",
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.01,
        #:mu_y                => 0.01,
        #:filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Topography parameters
        #---------------------------------------------------------------------------
        :lwarp                => true,
        :mount_type           => "real topography",
        :topo_database        => "./topography/ETOPO_2022_v1_30s_N90W180_bed.nc",
        :topo_geoid           => "./topography/ETOPO_2022_v1_30s_N90W180_geoid.nc",
        :read_topo_latmin     => 15.194166666,
        :read_topo_latmax     => 15.654166666,
        :read_topo_lonmin     => -61.5208333,
        :read_topo_lonmax     => -61.1208333,
        :read_topo_zone       => 20,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk", #"hdf5",
        :output_dir          => "./output/",
        :loverwrite_output   => true,
        :loutput_pert        => false,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
