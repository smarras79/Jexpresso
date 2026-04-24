function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 1D viscous Burgers equation (conservation form):
        #
        #     ∂q/∂t + ∂F(q)/∂x = ν ∂²q/∂x²,    F(q) = q²/2
        #
        # Time integration: IMEX Runge-Kutta ARS(2,3,2) from
        #   Ascher, Ruuth, Spiteri,
        #   "Implicit-explicit Runge-Kutta methods for time-dependent PDEs",
        #   Appl. Numer. Math. 25 (1997) 151-167.
        #
        #   - Advection (q²/2)_x  is treated EXPLICITLY
        #   - Diffusion ν q_xx    is treated IMPLICITLY (linear solve at each stage)
        #
        # The implicit split removes the parabolic Δt ∝ Δx² restriction and lets
        # the step size follow the (much milder) advective CFL condition.
        #---------------------------------------------------------------------------
        :ode_solver           => IMEX_ARS232(),
        :tend                 => 1.0,
        :Δt                   => 2.0e-3,
        :diagnostics_at_times => (0.1:0.1:1.0),
        :output_dir           => "./",
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 4,     # Polynomial order
        :lexact_integration  => false,
        :lsource             => false,
        :lperiodic_1d        => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #
        # :μ carries the kinematic viscosity ν that appears in the implicit
        # Laplacian built by imex_ars232_time_loop!. :lvisc is kept `true`
        # so that, if the user ever falls back to an explicit solver in the
        # same case, the existing viscous path still fires; the IMEX loop
        # temporarily toggles it off around each explicit-RHS evaluation.
        #---------------------------------------------------------------------------
        :lvisc               => true,
        :μ                   => [0.01],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => false,
        #---------------------------------------------------------------------------
        # Output formats
        #---------------------------------------------------------------------------
        :outformat           => "png",
        :loverwrite_output   => true,
        :output_dir          => "./output",
        #---------------------------------------------------------------------------
        # 1D (lread_gmsh => false): the grid is built by jexpresso
        #---------------------------------------------------------------------------
        :xmin                => 0.0,
        :xmax                => 1.0,
        :nelx                => 50,
    ) #Dict

    return inputs
end
