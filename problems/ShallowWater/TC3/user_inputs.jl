function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # TC3: Planetary Rossby Wave (Bishnu et al. 2024, Sec. 2.7)
        # Linearized rotating shallow water equations on a beta-plane
        # Domain [0, Lx] x [0, Ly] with sponge-layer absorbing BCs
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 10.0,
        :tinit                => 0.0,
        :tend                 => 1728000.0,   # 20 days in seconds
        :diagnostics_at_times => (0:86400:1728000),
        :case                 => "prw",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Filter for stabilization (CG SEM needs this for wave propagation)
        #---------------------------------------------------------------------------
        :lfilter              => true,
        :mu_x                 => 0.01,
        :mu_y                 => 0.01,
        :filter_type          => "erf",
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        # Generate mesh: gmsh -2 problems/ShallowWater/TC3/SWE_TC3.geo
        :gmsh_filename        => "./problems/ShallowWater/TC3/SWE_TC3.msh",
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
