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
        :tend                 => 100000.0,
        :diagnostics_at_times => (0:10000:100000),
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
    ) #Dict

    return inputs

end
