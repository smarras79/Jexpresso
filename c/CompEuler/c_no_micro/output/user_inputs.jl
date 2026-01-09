function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver            => CarpenterKennedy2N54(),
        :tend                 => 43200.0,
        #:ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :Δt                   => 0.01,
        :diagnostics_at_times => 0.0:600.0:43200.0,
        #:diagnostics_at_times => [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
        :case                 => "rtb",
        :lsource              => true,
        :lmoist               => false,    # Dry case (no moisture)
        :lprecip              => false,    # No microphysics
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations => [2, 3, 4, 5, 6],  # Momentum diffusion
        :μ                   => [0.0, 150.0, 150.0, 150.0, 150.0, 150.0], #horizontal viscosity constant for momentum
        :lfilter             => true,
        :mu_x                => 0.01,
        :mu_y                => 0.01,
        :filter_type         => "erf",
       # :sounding_file       => "./data_files/dry.dat",
      # :sounding_file       => "./data_files/strat.dat",
       #:sounding_file       => "./data_files/thermal.data",
       :sounding_file       => "./data_files/no_micro.data",
       # :sounding_file       => "./data_files/test_sounding.data",
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB20x20.msh", #for nop=4
       # :gmsh_filename       => "./meshes/gmsh_grids/coarse_4200_4.msh",
      #  :gmsh_filename       => "./meshes/gmsh_grids/squall2d.msh",
       :gmsh_filename       => "./meshes/gmsh_grids/c_4.msh",
       # :gmsh_filename       => "./meshes/gmsh_grids/square_UNSTR_20el.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk", #"hdf5",
        :loverwrite_output   => true,
        :output_dir          => "./c",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
