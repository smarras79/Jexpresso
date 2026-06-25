function user_inputs()

    #--------------------------------------------------------------------------
    # IMEX Runge-Kutta ARS(2,3,2) Butcher tableaux
    #
    #   Ascher, Ruuth, Spiteri, "Implicit-explicit Runge-Kutta methods for
    #   time-dependent PDEs", Appl. Numer. Math. 25 (1997) 151-167.
    #
    # (A_RK,b_RK,c_RK)             -> explicit part (advection)
    # (A_RK_tilde,b_RK_tilde,c_RK_tilde) -> implicit DIRK part (diffusion)
    #--------------------------------------------------------------------------
    a_32 = 1.0 / 6.0 * (3.0 + 2.0 * sqrt(2.0))

    A_RK = zeros(TFloat, 3, 3)
    A_RK[2, 1] = TFloat(2.0 - sqrt(2.0))
    A_RK[3, 1] = TFloat(1.0 - a_32)
    A_RK[3, 2] = TFloat(a_32)

    b_RK = TFloat[ 1.0 / (2.0 * sqrt(2.0)),
                   1.0 / (2.0 * sqrt(2.0)),
                   1.0 - 1.0 / sqrt(2.0) ]

    c_RK = TFloat[ 0.0, 2.0 - sqrt(2.0), 1.0 ]

    A_RK_tilde = zeros(TFloat, 3, 3)
    A_RK_tilde[2, 1] = TFloat(1.0 - 1.0 / sqrt(2.0))
    A_RK_tilde[2, 2] = TFloat(1.0 - 1.0 / sqrt(2.0))
    A_RK_tilde[3, 1] = TFloat(1.0 / (2.0 * sqrt(2.0)))
    A_RK_tilde[3, 2] = TFloat(1.0 / (2.0 * sqrt(2.0)))
    A_RK_tilde[3, 3] = TFloat(1.0 - 1.0 / sqrt(2.0))

    b_RK_tilde = TFloat[ 1.0 / (2.0 * sqrt(2.0)),
                         1.0 / (2.0 * sqrt(2.0)),
                         1.0 - 1.0 / sqrt(2.0) ]

    c_RK_tilde = TFloat[ 0.0, 2.0 - sqrt(2.0), 1.0 ]

    #--------------------------------------------------------------------------
    # Explicit source S(u): the full explicit RHS, evaluated through the
    # standard Jexpresso `rhs!` (advection plus the explicit viscous term that
    # `:lvisc => true` enables). The implicit operator L below provides the
    # additional, unconditionally-stable diffusion used inside the per-stage
    # solves.
    #--------------------------------------------------------------------------
    function S_fun!(s_j, u, time, params, sem)
        rhs!(s_j, u, params, time)
        s_j .= params.RHS
    end

    #--------------------------------------------------------------------------
    # Implicit operator  L = -μ M⁻¹ K   (K = positive 2D Galerkin stiffness).
    #
    # It is linear and time-independent, so the sparse assembly
    # (build_laplace_matrix + DSS_laplace_sparse, the native Jexpresso
    # element-matrix infrastructure) only has to run once. The assembled
    # operator is cached in the closure so build_L / L_fun! never re-assemble.
    #--------------------------------------------------------------------------
    L_cache = Ref{Any}(nothing)

    function _imex_L(params)
        if L_cache[] === nothing
            SD      = params.SD
            basis   = params.basis
            ω       = params.ω
            mesh    = params.mesh
            metrics = params.metrics
            μ       = params.inputs[:μ]
            N       = params.inputs[:nop]
            Q       = N

            Le       = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh.nelem,
                                            mesh, metrics, N, Q, TFloat)
            L_global = DSS_laplace_sparse(mesh, Le)
            Minv     = params.Minv
            # NSD_2D build_laplace_matrix returns the positive stiffness K, so
            # the diffusion operator ν Δ ↦ -μ M⁻¹ K picks up a minus sign.
            L_cache[] = -μ[1] * (Minv .* L_global)
        end
        return L_cache[]
    end

    # Apply the fast (diffusive) operator: l_j = L u, via cached sparse mul!.
    function L_fun!(l_j, u, time, params)
        mul!(l_j, _imex_L(params), u)
    end

    # Assemble (and cache) the implicit operator L = -μ M⁻¹ K.
    function build_L(u, time, params)
        return _imex_L(params)
    end

    # Per-stage linear solve. The implicit operator is a well-conditioned
    # sparse SPD-like Helmholtz operator I - λ L, so a direct sparse solve is
    # used here.
    function lsolve(L_curr, b)
        return L_curr \ b
    end

    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D Burgers Riemann problem (flux form), Sec. 4.1 of
        #
        #     ∂ₜ u + ∇·(½ u² v) = 0,    v = (1, 1)
        #
        # IMEX (implicit-explicit) time integration of the viscous-regularized
        # Burgers equation:
        #
        #   - Advection ∇·(½ u² v)   is treated EXPLICITLY via rhs!  (S_fun!)
        #   - Diffusion ν Δu         is treated IMPLICITLY via the sparse
        #                            operator L = -μ M⁻¹ K            (L_fun!)
        #
        # Time scheme: ARS(2,3,2) additive Runge-Kutta (tableaux above).
        #
        # This is the "_sl" companion of Burgers/case2d: it exercises
        # Jexpresso's native IMEX integrator (`:ode_solver => IMEX()`,
        # kernel/solvers/IMEXTimeIntegrators.jl) rather than the explicit
        # OrdinaryDiffEq path. The physics, mesh and initial condition are
        # identical to case2d.
        #
        # Mesh: a doubly-periodic quadrilateral gmsh mesh. The Laplacian
        # stiffness inherits the periodicity from mesh.connijk, so no extra BC
        # wiring is needed.
        #---------------------------------------------------------------------------
        :ode_solver           => IMEX(),
        :tinit                => 0.0,
        :tend                 => 0.5,
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => (0.1:0.1:0.5),
        :output_dir           => "./output",
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
        # :μ carries the kinematic viscosity ν. It is used both by the explicit
        # viscous RHS (enabled by :lvisc => true) and to assemble the implicit
        # Laplacian L = -μ M⁻¹ K used to stabilize the per-stage solves. Set
        # :μ => [0.0, 0.0] to recover the pure inviscid problem.
        #---------------------------------------------------------------------------
        :lvisc               => true,
        :μ                   => [1.0e-2, 1.0e-2],
        #---------------------------------------------------------------------------
        # Mesh parameters and files (doubly-periodic quadrilateral mesh)
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh",
        #---------------------------------------------------------------------------
        # Output formats
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        #---------------------------------------------------------------------------
        # Matrix assembly: build the sparse Galerkin Laplacian used by the
        # implicit operator.
        #---------------------------------------------------------------------------
        :ldss_laplace        => true,
        :lsparse             => true,
        :ldss_differentiation => false,
        #---------------------------------------------------------------------------
        # IMEX method configuration (consumed by imex_time_loop!)
        #---------------------------------------------------------------------------
        :method              => "RK",
        :delta               => 1,      # 0 -> explicit, 1 -> IMEX
        :k                   => 3,      # # of RK stages
        :coeff               => Dict(
                                    :A_RK       => A_RK,
                                    :b_RK       => b_RK,
                                    :c_RK       => c_RK,
                                    :A_RK_tilde => A_RK_tilde,
                                    :b_RK_tilde => b_RK_tilde,
                                    :c_RK_tilde => c_RK_tilde,
                                ),
        :S_fun               => S_fun!,
        :L_fun               => L_fun!,
        :build_L             => build_L,
        :upd_L               => false,
        :lsolve              => lsolve,
        :solver_precision    => TFloat,
        :nl_precision        => TFloat,
    ) #Dict

    return inputs
end
