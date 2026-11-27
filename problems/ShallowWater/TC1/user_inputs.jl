function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # TC1: Coastal Kelvin Wave (Bishnu et al. 2024, Sec. 2.5)
        # Linearized rotating shallow water equations on an f-plane
        # Coast (solid wall) at x=0, open boundary at x=Lx, periodic in y
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 50.0,    # paper DGSEM dt (Bishnu et al. 2024, Sec. 6)
        :tinit                => 0.0,
        :tend                 => 400000.0,           # ~8 periods of the longest mode (Ly/c = 5e4 s)
        :diagnostics_at_times => (0:50000:400000),   # sparse VTK, for a quick look only
        :case                 => "kelvin",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 10,      # paper: N = 5 elements, P = 10
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Filter: off. With it on (erf, mu = 0.05) it was ~99% of TC1's error
        # (~1e-7 vs ~1e-9 without, at 3-9 h) and not needed for stability; the
        # paper's runs are unfiltered too.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        :mu_x                 => 0.05,
        :mu_y                 => 0.05,
        :filter_type          => "erf",
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        # Mesh: meshes/gmsh_grids/SWE_TC1_coastal.geo (gmsh -2 -format msh41 SWE_TC1_coastal.geo -o SWE_TC1_coastal_5x5.msh)
        :gmsh_filename        => "./meshes/gmsh_grids/SWE_TC1_coastal_5x5.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
        # ML snapshot recording (RF-QRC proposal, Task 1.1), as in
        # SWsphere_ScottPolvani_POD: the POD recorder samples η, u, v at exactly
        # uniform times (the integrator stops on them) and writes the raw
        # snapshots; node coordinates come from the pod_<field> VTK of the same
        # run (same nodes, same order). The decomposition itself is a by-product.
        #---------------------------------------------------------------------------
        :lpod                 => true,
        :pod_fields           => ["η", "u", "v"],
        :pod_tstart           => 0.0,
        :pod_tend             => 400000.0,
        :pod_nsnapshots       => 800,      # intervals -> 801 snapshots, Δt_save = 500 s
        :pod_nmodes           => 0,        # keep every mode (the snapshots are small)
        :pod_write_snapshots  => true,     # the ML data
        :pod_write_vtk        => true,     # carries the node coordinates for the snapshots
        :pod_write_png        => false,
        :pod_write_data       => true,
        :pod_time_scale       => 1.0/3600.0,
        :pod_time_label       => "t [h]",
    ) #Dict

    return inputs

end
