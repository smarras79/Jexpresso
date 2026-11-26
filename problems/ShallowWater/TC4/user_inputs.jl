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
    ) #Dict

    return inputs

end
