function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 1D inviscid Burgers equation, stabilized with the Marras et al.
        # Dynamic SGS (DSGS) artificial viscosity, as CompEuler/sod1d:
        #
        #     ∂q/∂t + ∂(q²/2)/∂x = ∂ₓ(ν_DSGS ∂ₓq),    q(x,0) = sin(2πx)
        #
        # on the periodic unit interval. A stationary shock forms at x = 0.5
        # at t = 1/(2π) ≈ 0.159; the exact entropy solution is overlaid on
        # every frame (user_analytic.jl). Without DSGS the spectral element
        # solution Gibbs-oscillates at the shock (see Burgers/case1 for the
        # constant-viscosity version).
        #
        # The coefficient is the scalar branch of the 1D DSGS kernel
        # (_dsgs_1d_burgers! in src/kernel/physics/SGS.jl): the element
        # residual of the one equation, capped by Cmax·Δ·max|q|, |q| being
        # the Burgers characteristic speed.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK53(),
        :tend                 => 0.5,
        :Δt                   => 2.0e-4,
        :diagnostics_at_times => (0:0.05:0.5),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => false,
        :lperiodic_1d         => true,
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => DSGS(),
        :dsgs_sensor          => "residual",  # element-wise strong residual (DSGS.md §1.2)
        # Multiplier on the DSGS coefficient: 1.0 = Marras value.
        :μ                    => [1.0],
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
