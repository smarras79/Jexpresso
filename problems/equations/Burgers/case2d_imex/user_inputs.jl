function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D viscous Burgers equation (conservation form):
        #
        #     ∂q/∂t + ∂F/∂x + ∂G/∂y = ν (∂²q/∂x² + ∂²q/∂y²),
        #         F(q) = q²/2,  G(q) = q²/2
        #
        # Time integration: IMEX Runge-Kutta ARS(2,3,2) from
        #   Ascher, Ruuth, Spiteri,
        #   "Implicit-explicit Runge-Kutta methods for time-dependent PDEs",
        #   Appl. Numer. Math. 25 (1997) 151-167.
        #
        #   - Advection (q²/2)_x + (q²/2)_y  is treated EXPLICITLY via rhs!
        #   - Diffusion ν Δq                 is treated IMPLICITLY (sparse solve)
        #
        # The implicit split removes the parabolic Δt ∝ h² restriction and lets
        # the step size follow the (much milder) advective CFL condition.
        #
        # Mesh: a doubly-periodic quadrilateral gmsh mesh. The Laplacian
        # stiffness assembled by imex_ars232_time_loop! inherits the
        # periodicity from mesh.connijk, so no extra BC wiring is needed.
        #---------------------------------------------------------------------------
        :ode_solver           => IMEX_ARS232(),
        :tend                 => 0.5,
        :Δt                   => 2.0e-3,
        :diagnostics_at_times => (0.05:0.05:0.5),
        :output_dir           => "./",
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 4,     # Polynomial order
        :lexact_integration  => false,
        :lsource             => false,
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
        :ivisc_equations     => [1],
        :μ                   => [0.01],
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #   Reuse the doubly-periodic quadrilateral mesh shipped with the
        #   2D AdvDiff demo. Drop in any other periodic 2D gmsh mesh as long
        #   as the polynomial order matches :nop above.
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10_periodic.msh",
        #---------------------------------------------------------------------------
        # Output formats
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        :output_dir          => "./output",
    ) #Dict

    return inputs
end
