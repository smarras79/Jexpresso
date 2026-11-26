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
        :tend                 => 100000.0,
        :diagnostics_at_times => (0:10000:100000),
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
    ) #Dict
    return inputs
end
