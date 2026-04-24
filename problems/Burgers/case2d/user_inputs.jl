function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D Burgers Riemann problem (flux form), Sec. 4.1 of
        #
        #     ∂ₜ u + ∇·(½ u² v) = 0,    v = (1, 1)
        #
        # This is the inviscid test; the implicit side of the IMEX split
        # covers a (possibly zero) viscous regularization ν Δu. A small
        # ν > 0 is recommended because a pure CG-SEM discretization does
        # not limit Gibbs oscillations at the initial discontinuities —
        # :μ below acts as vanishing viscosity. Set :μ => [0.0] to
        # reproduce the purely inviscid problem from the paper.
        #
        # Time integration: IMEX Runge-Kutta ARS(2,3,2) from
        #   Ascher, Ruuth, Spiteri,
        #   "Implicit-explicit Runge-Kutta methods for time-dependent PDEs",
        #   Appl. Numer. Math. 25 (1997) 151-167.
        #
        #   - Advection ∇·(½ u² v)   is treated EXPLICITLY via rhs!
        #   - Diffusion ν Δu         is treated IMPLICITLY (sparse solve)
        #
        # Mesh: a doubly-periodic quadrilateral gmsh mesh. The Laplacian
        # stiffness assembled by imex_ars232_time_loop! inherits the
        # periodicity from mesh.connijk, so no extra BC wiring is needed.
        # Initial discontinuities sit at the domain midpoints, i.e. the
        # reference paper's x = y = 0.5 lines when the mesh is the unit
        # square.
        #---------------------------------------------------------------------------
        #:ode_solver           => IMEX_ARS232(),
        :ode_solver           => CarpenterKennedy2N54(),
        :tend                 => 0.5,
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => (0.05:0.05:0.5),
        :output_dir           => "./",
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl", # Choice: "lgl", "cg", "cgl"
        :nop                 => 7,     # Polynomial order
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
        :μ                   => [1.0e-2, 1.0e-2],   # set to 0.0 for the pure inviscid case in the paper
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #   Reuse the doubly-periodic quadrilateral mesh shipped with the
        #   2D AdvDiff demo. Drop in any other periodic 2D gmsh mesh as long
        #   as the polynomial order matches :nop above.
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh",
        #---------------------------------------------------------------------------
        # Output formats
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        :output_dir          => "./output",
    ) #Dict

    return inputs
end
