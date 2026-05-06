function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        #:Δt                   => 0.6, for 1600
        #:Δt                   => 0.4, for 800
        :Δt                   => 0.3, #for 400
        :tinit                => 0.0,
        :tend                 => 28800.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        #:restart_input_file_path => "./output/CompEuler/theta/output-19Nov2023-115126",
        :diagnostics_at_times => (0:60:28800),
        :case                 => "rtb",
        :lsource              => true, 
        #:lmoist               => true,
        #:lprecip              => true,
        :lsponge              => true,
        :zsponge              => 20000.0,
        :SOL_VARS_TYPE        => PERT(),
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :nop                 => 4,      # Polynomial order
       # :nop_laguerre        => 24,
          #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc       => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => (2,3,4),
        :μ => [0.0, 25.0, 25.0, 0.0],
        #:visc_model  => SMAG()
        #:μ           => [0.0, 20.0, 20.0, 30.0, 30.0, 30.0], #horizontal viscosity constant for momentum
       # :μ           => [0.0, 400.0, 400.0, 400.0], #horizontal viscosity constant for momentum
      # :μ           => [0.0, 200.0, 200.0, 200.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/agnesi240kmX30km_coarse.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10_laguerre_top.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/agnesi-120kmx30km-hm5000.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_RTB.msh",
       # :gmsh_filename         => "./meshes/gmsh_grids/m400_4.msh",
       :gmsh_filename         => "./meshes/gmsh_grids/snwp.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/agnesi240kmX30km_coarse_laguerreTopLateral.msh",
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        :lwarp               => true,
        :mount_type          => "agnesi",
        :a_mount             => 10000.0,
        :h_mount  => 1.0,   # currently 100.0 — needs to be 2.5 km
        :c_mount  => 40000.0,  # currently 0.0 — mountain center at 40 km   
        #---------------------------------------------------------------------------
        # Soundings and data files
        #---------------------------------------------------------------------------
        #:sounding_file       => "./data_files/test_hill_hydro.data", #hydrostatic balanced test_hill_data.
        :sounding_file       => "./data_files/SNWP.data", #no_shear data.
        #:sounding_file       => "./data_files/pap.data", #main sounding data
       #  :sounding_file       => "./data_files/test_hill_hydro_noqv.data", 
       # :sounding_file       => "./data_files/input_sounding", #test_hill.data minus qv
        #:sounding_file       => "./data_files/sounding-SAM-new.dat", 
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => true,
        :mu_x                => 0.25,
        :mu_y                => 0.25,
        :filter_type         => "erf",  ##default is erf, use either "erf" for Boyd-Vandeven,"exp" for Warburton Exponential filter, or "quad" for Fischer quadratic filter
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        #:output_dir          => "./truly_dry",
        :output_dir          => "./SNWP",
        :loverwrite_output   => true,
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
