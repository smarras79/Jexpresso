# `user_inputs.jl`
Every user-defined problem has its own namelist/input file called `user_inputs.jl`.
A sample user_inputs.jl file with all possible entries (notice that not all entries are necessary).
Use `#` to comment a line. When a line is commented or not explicitaly give, default values will be used. 
The entries (with default values) are defined in `src/io/mod_inputs.jl`.


```
function user_inputs()
    
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        #:Δt                   => 0.02,
        :Δt                   => 0.4,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        :restart_input_file_path => "./output/CompEuler/theta/",
        :diagnostics_at_times => (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
        :lsource              => true, 
        #:backend              => MetalBackend(),
        #:SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true,                   #false by default
        :ivisc_equations      => [1, 2, 3, 4],           # use [], not ()!
        :μ                   => [0.0, 20.0, 20.0, 60.0], # use [], not ()!
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB20x20.msh",	
        #---------------------------------------------------------------------------
	#
	# Topography parameters:
	# Real topography can be extreacted from the ETOPO data by defininging
	# the domain of interest by means of the lat-lon extrema.
	#
	# The user should download the ETOPO global relief model by NOAA:
	#
	# https://www.ncei.noaa.gov/products/etopo-global-relief-model
	#
	# You need the geoid file *_geoid.nc
	# 
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
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.01,
        #:mu_y                => 0.01,
        #:filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        :lwrite_initial      => false,
        :output_dir          => "./output",
        #:output_dir          => "./test/CI-run",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
```