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
        visc_coeff = params.visc_coeff
        N          = nop
        Q          = N
        backend    = CPU()

        Le = KernelAbstractions.zeros(backend, TFloat, 1, 1)
        L  = KernelAbstractions.zeros(backend, TFloat, 1, 1)

        Le = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh, metrics,
                                  N, Q, TFloat)

        L = DSS_laplace_sparse(mesh, Le)
        L = - visc_coeff[1] * L
        
        return L
    end

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
#        :ode_solver           => IMEX_ARS232(),
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
        :delta              => 1,
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
