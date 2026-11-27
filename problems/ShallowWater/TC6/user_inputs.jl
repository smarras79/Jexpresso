function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # TC6: Non-linear manufactured solution (Bishnu et al. 2024, Sec. 2.10)
        # Non-linear rotating SWE, conservative form, f-plane, doubly periodic.
        # Paper DGSEM resolution (Sec. 6): N = 5 elements, P = 10, dt = 30 s.
        # Same 1e7 m doubly periodic domain as TC2, so it reuses TC2's mesh.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 30.0,     # paper DGSEM dt
        :tinit                => 0.0,
        :tend                 => 566400.0,           # ~8 periods (T = 19.6 h)
        :diagnostics_at_times => (0:141600:566400),  # sparse VTK, for a quick look only
        :case                 => "mms",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 10,       # paper: N = 5 elements, P = 10
        #---------------------------------------------------------------------------
        # Physical parameters/constants: inviscid (paper Eq. 39/41)
        #---------------------------------------------------------------------------
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
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
        # ML snapshot recording (RF-QRC proposal, Task 1.1), as in TC1/TC2: the POD
        # recorder samples η, u, v at exactly uniform times; node coordinates
        # come from the pod_<field> VTK of the same run.
        #---------------------------------------------------------------------------
        :lpod                 => true,
        :pod_fields           => ["η", "u", "v"],
        :pod_tstart           => 0.0,
        :pod_tend             => 566400.0,
        :pod_nsnapshots       => 944,      # intervals -> 945 snapshots, Δt_save = 600 s
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
