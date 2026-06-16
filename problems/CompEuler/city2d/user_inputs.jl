function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 0.004,
        :tinit                => 0.0,
        :tend                 => 600.0,
        :diagnostics_at_times => (0:10:600),
        :restart_time         => 50000,
        :lrestart             => false,
        :case                 => "city2d",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Sponge: weak top sponge of thickness 500 m (ymax=1000 -> zsponge=500)
        #---------------------------------------------------------------------------
        :lsponge              => true,
        :zsponge              => 500.0,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 4,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc           => true,
        :visc_model      => SMAG(),
        :energy_equation => "theta",
        :μ               => [0.0, 5.0, 5.0, 5.0],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh    => true,
        #:gmsh_filename => "./meshes/gmsh_grids/city2d.msh",
        :gmsh_filename => "./meshes/gmsh_grids/city2d_transfinite.msh",
        #:gmsh_filename => "./meshes/gmsh_grids/city2d_TFI.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :lstep_heartbeat   => true, #comment it if you don't want the steps printing
        :outformat         => "vtk",
        :loverwrite_output => true,
        :lwrite_initial    => true,
        :output_dir        => "./output",
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine => false,
        :init_refine_lvl => 1,
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :ladapt    => false,
        :amr_freq  => 200,
        :amr_max_level => 2,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs

end
