function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # TC2: Inertia-Gravity Wave (Bishnu et al. 2024, Sec. 2.3)
        # Linearized rotating shallow water equations on an f-plane
        # Doubly periodic domain [0, Lx] x [0, Ly]
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 25.0,    # paper DGSEM dt (Bishnu et al. 2024, Sec. 6)
        :tinit                => 0.0,
        :tend                 => 266400.0,           # ~74 h, as in the proposal's preliminary TC2
        :diagnostics_at_times => (0:66600:266400),   # sparse VTK, for a quick look only
        :case                 => "igw",
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
        # Filter for stabilization (CG SEM needs this for wave propagation)
        #---------------------------------------------------------------------------
#        :lfilter              => true,
        :mu_x                 => 0.05,
        :mu_y                 => 0.05,
        :filter_type          => "erf",
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        # Generate mesh: gmsh -2 problems/ShallowWater/TC2/SWE_TC2_periodic.geo
        :gmsh_filename        => "./meshes/gmsh_grids/SWE_TC2_periodic_5x5.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
        # ML snapshot recording (RF-QRC proposal, Task 1.1), as in TC1: the POD
        # recorder samples η, u, v at exactly uniform times; node coordinates
        # come from the pod_<field> VTK of the same run.
        #---------------------------------------------------------------------------
        :lpod                 => true,
        :pod_fields           => ["η", "u", "v"],
        :pod_tstart           => 0.0,
        :pod_tend             => 266400.0,
        :pod_nsnapshots       => 592,      # intervals -> 593 snapshots, Δt_save = 450 s
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
