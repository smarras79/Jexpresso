function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 1D CompEuler: Sod's shock-tube Riemann problem, stabilized with the
        # Marras et al. Dynamic SGS (DSGS) artificial viscosity.
        #
        # Initial condition (piecewise constant about x = 0.5):
        #     (ρ, u, p) = (1.000, 0, 1.0)   for x < 0.5
        #     (ρ, u, p) = (0.125, 0, 0.1)   for x ≥ 0.5
        #
        # The reference solution at t ≈ 0.2 shows a left rarefaction, a contact
        # discontinuity, and a right-moving shock. DSGS is required because the
        # raw spectral element discretisation Gibbs-oscillates around the shock.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK53(),
        :tend                 => 0.2,
        :Δt                   => 1.0e-4,
        :diagnostics_at_times => (0:0.02:0.2),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => false,
        :lperiodic_1d         => false,
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => DSGS(),
        # Per-equation multiplier on the DSGS coefficient: 1.0 = Marras
        # value, 0.0 turns it off on that equation, in (0, 1) throttles.
        # Mass diffusion is kept on for 1D shock stabilization (Marras's
        # 1D analysis includes it; switching it off lets the Sod shock
        # ring).
        :μ                    => [1.0, 1.0, 1.0],
        #---------------------------------------------------------------------------
        # Mesh parameters and files
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => 0.0,
        :xmax                 => 1.0,
        :nelx                 => 100,
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs

end
