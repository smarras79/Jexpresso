function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # TC4: Topographic Rossby Wave (Bishnu et al. 2024, Sec. 2.8)
        # Linearized rotating SWE on an f-plane with a sloping bottom D(y) = H0 + a0*y.
        # Same domain, mesh, resolution, IC and wall BCs as TC3 (paper Sec. 2.8).
        # Domain [0, Lx] x [0, Ly] with sponge-layer absorbing BCs
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 10.0,
        :tinit                => 0.0,
        :tend                 => 1728000.0,   # 20 days in seconds
        :diagnostics_at_times => (0:86400:1728000),
        :case                 => "trw",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 7,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Filter: off (as TC1, TC2, TC6 and the paper). Was on (erf, mu = 0.01).
        #---------------------------------------------------------------------------
        :lfilter              => false,
        :mu_x                 => 0.01,
        :mu_y                 => 0.01,
        :filter_type          => "erf",
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        # Same mesh as TC3
        :gmsh_filename        => "./meshes/gmsh_grids/SWE_TC3_10x10.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
        #---------------------------------------------------------------------------
        # ML snapshot recording (RF-QRC proposal, Task 1.1), as in TC1/TC2: the POD
        # recorder samples η, u, v at exactly uniform times; node coordinates
        # come from the pod_<field> VTK of the same run.
        #---------------------------------------------------------------------------
        :lpod                 => true,
        :pod_fields           => ["η", "u", "v"],
        :pod_tstart           => 0.0,
        :pod_tend             => 1727000.0,
        :pod_nsnapshots       => 1570,      # intervals -> 1571 snapshots, Δt_save = 1100 s (window ends 1000 s before :tend)
        :pod_nmodes           => 0,
        :pod_write_snapshots  => true,     # the ML data
        :pod_write_vtk        => true,     # carries the node coordinates for the snapshots
        :pod_write_png        => false,
        :pod_write_data       => true,
        :pod_time_scale       => 1.0/3600.0,
        :pod_time_label       => "t [h]",
    ) #Dict

    return inputs

end
