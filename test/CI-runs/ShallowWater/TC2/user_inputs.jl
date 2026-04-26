function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # CI version: reduced tend for fast execution (2 time steps)
        #
        # TC2: Inertia-Gravity Wave
        # Linearized rotating shallow water equations on an f-plane
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 50.0,
        :tinit                => 0.0,
        :tend                 => 100.0,
        :diagnostics_at_times => (100.0,),
        :case                 => "igw",
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
        # Filter for stabilization
        #---------------------------------------------------------------------------
        :mu_x                 => 0.05,
        :mu_y                 => 0.05,
        :filter_type          => "erf",
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/SWE_TC2_periodic.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat            => "hdf5",
        :loverwrite_output    => true,
        :lwrite_initial       => false,
        :output_dir           => "none",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
