function user_inputs()
    # Runge-Kutta
    a_32 = 1. / 6. * (3. + 2. * sqrt(2.))
    A_RK = zeros(3, 3)
    A_RK[1, 1] = 0.
    A_RK[1, 2] = 0.
    A_RK[1, 3] = 0.
    A_RK[2, 1] = 2. - sqrt(2.)
    A_RK[2, 2] = 0.
    A_RK[2, 3] = 0.
    A_RK[3, 1] = 1. - a_32
    A_RK[3, 2] = a_32
    A_RK[3, 3] = 0.

    b_RK = Array{Float64, 1}(undef, 3)
    b_RK[1] = 1. / (2. * sqrt(2.))
    b_RK[2] = 1. / (2. * sqrt(2.))
    b_RK[3] = 1. - 1. / sqrt(2.)

    c_RK = Array{Float64, 1}(undef, 3)
    c_RK[1] = 0.
    c_RK[2] = 2. - sqrt(2.)
    c_RK[3] = 1.

    A_RK_tilde = zeros(3, 3)
    A_RK_tilde[1, 1] = 0.
    A_RK_tilde[1, 2] = 0.
    A_RK_tilde[1, 3] = 0.
    A_RK_tilde[2, 1] = 1. - 1. / sqrt(2.)
    A_RK_tilde[2, 2] = 1. - 1. / sqrt(2.)
    A_RK_tilde[2, 3] = 0.
    A_RK_tilde[3, 1] = 1. / (2. * sqrt(2.))
    A_RK_tilde[3, 2] = 1. / (2. * sqrt(2.))
    A_RK_tilde[3, 3] = 1. - 1. / sqrt(2.)

    b_RK_tilde = Array{Float64, 1}(undef, 3)
    b_RK_tilde[1] = 1. / (2. * sqrt(2.))
    b_RK_tilde[2] = 1. / (2. * sqrt(2.))
    b_RK_tilde[3] = 1. - 1. / sqrt(2.)

    c_RK_tilde = Array{Float64, 1}(undef, 3)
    c_RK_tilde[1] = 0.
    c_RK_tilde[2] = 2. - sqrt(2.)
    c_RK_tilde[3] = 1.

    # Polynomial order
    nop = 4

    # Solver parameters
    solver_par = Dict(:restart  => true,
                      :memory   => 10,
                      :verbose  => 1,
                      :atol     => 1.e-10,
                      :rtol     => 1.e-10,
                      :itmax    => 100,
                      :prec     => SmoothedAggregationPreconBuilder()
                      )

    # Preconditioner parameters
    prec_sp = Dict(
        :maxiter      => 1,
        :abstol       => 1e-8,
        :precision    => Float32,
        :prec_type    => "AMG",
        )

    # Source function
    function S_fun!(s_j, u, time, params, sem)
        rhs!(s_j, u, params, time)

        s_j .= params.RHS
    end

    # Fast waves operator
    function L_fun!(l_j, u, time, params)
        SD    = params.SD
        QT    = params.QT
        AD    = params.AD
        neqs  = params.neqs
        ngl   = params.mesh.ngl
        nelem = params.mesh.nelem

        resetRHSToZero_viscous!(params, SD)
        
        viscous_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, SD)
        
        # @info "start DSS_rhs_viscous"
        if inputs[:ladapt] == true
            DSS_nc_gather_rhs!(params.RHS_visc, SD, QT, params.rhs_diff_el, params.mesh.connijk, params.mesh.poin_in_edge, params.mesh.non_conforming_facets,
                               params.mesh.non_conforming_facets_parents_ghost, params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost, params.mesh.pgip_owner, ngl-1, neqs, params.interp)
        end
        DSS_rhs!(params.RHS_visc, params.rhs_diff_el, params.mesh.connijk, nelem, ngl, neqs, SD, AD)

        DSS_global_RHS!(@view(params.RHS_visc[:,:]), params.pM, params.neqs)

        l_j .= params.RHS_visc
    end

    # Bcs application
    function bcs_fun!(u, L, time, params, sem, qp)
        apply_boundary_conditions_lin_solve!(L, time, params.qp.qe,
                                             params.mesh.x, params.mesh.y, params.mesh.z,
                                             params.metrics.nx,
                                             params.metrics.ny,
                                             params.metrics.nz,
                                             sem.mesh.npoin, params.mesh.npoin_linear, 
                                             params.mesh.poin_in_bdy_edge,
                                             params.mesh.poin_in_bdy_face,
                                             params.mesh.nedges_bdy,
                                             params.mesh.nfaces_bdy,
                                             params.mesh.ngl, params.mesh.ngr,
                                             params.mesh.nelem_semi_inf,
                                             params.basis.ψ, params.basis.dψ,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             u, 0.0, params.ubdy,
                                             params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem,
                                             params.mesh.bdy_edge_type,
                                             params.ω, qp.neqs, params.inputs, params.AD, sem.mesh.SD)
    end

    # Building fast waves operator
    function build_L(u, time, params)
        SD         = params.SD
        basis      = params.basis
        ω          = params.ω
        mesh       = params.mesh
        metrics    = params.metrics
        μ          = 1.0e-2
        N          = nop
        Q          = N
        backend    = CPU()

        Le = KernelAbstractions.zeros(backend, TFloat, 1, 1)
        L  = KernelAbstractions.zeros(backend, TFloat, 1, 1)

        Le = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh.nelem, mesh, metrics,
                                  N, Q, TFloat)

        L = DSS_laplace_sparse(mesh, Le)
        L = - μ * L
        
        return L
    end

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
#        :ode_solver           => IMEX_ARS232(),
        :tend                 => 0.5,
        :Δt                   => 1.0e-4,
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
        #---------------------------------------------------------------------------
        #Building matrices
        #---------------------------------------------------------------------------
        :ldss_laplace        => true,
        :lsparse             => true,
        :ldss_differentiation => false,#true, # with true I get error (DSS_generic_matrix not defined)
        #---------------------------------------------------------------------------
        # IMEX method
        #---------------------------------------------------------------------------
        :method             => "RK",
        :delta              => 0,
        :k                  => 3,
        :coeff              => Dict(
                                   # IMEX RK
                                   :A_RK        => A_RK,
                                   :b_RK        => b_RK,
                                   :c_RK        => c_RK,
                                   :A_RK_tilde  => A_RK_tilde,
                                   :b_RK_tilde  => b_RK_tilde,
                                   :c_RK_tilde  => c_RK_tilde,
                               ),
        :lsolver            => nothing,#"GMRES",#LinearSolve.KrylovJL_GMRES(),
        :sp                 => solver_par,
        :prec_sp            => prec_sp,
        :S_fun              => S_fun!,
        :L_fun              => L_fun!,
        :bcs_fun            => bcs_fun!,
        :upd_L              => false,
        :build_L            => build_L,
        #---------------------------------------------------------------------------
        # Matrix storage for the IMEX implicit operator
        #   :matrix_free  - (default) rebuild L and the preconditioner on demand
        #   :assembled    - store L in optimal sparse (CSC) storage together
        #                   with mixed-precision copies and a cached
        #                   preconditioner. Required to use ILU / tuned AMG
        #                   preconditioners or mixed-precision strategies.
        #---------------------------------------------------------------------------
        :matrix_storage     => :assembled,
    ) #Dict

    return inputs
end
