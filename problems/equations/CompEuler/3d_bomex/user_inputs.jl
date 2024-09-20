function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.1,
        :tinit                => 0.0,
        :tend                 => 10000,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        :restart_input_file_path => "./output/CompEuler/theta/output-19Nov2023-115126",
        :diagnostics_at_times => (100:100:10000),
        # :ndiagnostics_outputs => 11,
        :case                 => "rtb",
        :lbomex               => true,
        :lsource              => true, 
        :backend              => CUDABackend(),
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :visc_model           => "smag",
        :ivisc_equations      => [1, 2, 3, 4, 5, 6, 7],
        # smagorinsky, cs = 0.23, input cs^2 for momentum cs^2/Pr for other equations, where Pr = 1/3
        :μ                    => [0.1587, 0.0529, 0.0529, 0.0529, 0.1587, 0.1587, 0.1587], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh           => true, #If false, a 1D problem will be enforced
        # :gmsh_filename        => "./meshes/gmsh_grids/hexa_BOMEX-5x1x19.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_BOMEX-16x16x19.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_BOMEX-32x32x38.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x1x10.msh",
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
        :loutput_pert        => false,  #this is only implemented for VTK for now
        :lwrite_initial      => true,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
