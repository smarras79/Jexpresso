function user_inputs()
    # Coefficients of the method
    # multistep
    alpha = Array{Float64, 1}(undef, 2)
    alpha[1] = 4. / 3.
    alpha[2] = - 1. / 3.
    beta = Array{Float64, 1}(undef, 2)
    beta[1] = 2.
    beta[2] = - 1.

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
        )

    # Source function
    function S_fun!(s_j, u, time, params)
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

        Le = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh.nelem, mesh, metrics,
                                  N, Q, TFloat)

        L = DSS_laplace_sparse(mesh, Le)
        assemble_diffusion_matrix_threaded!(mesh, Le)
        L = - visc_coeff[1] * L

        return L
    end

    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 1.0, #2π,
        :Δt                   => 0.001,#8.75e-4,
        :Δt_expl              => 0.0001,#8.75e-4,
        :ode_solver           => SSPRK54(),
        :diagnostics_at_times => (4.0),
        :output_dir          => "./output/",
#        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :nop                 => nop,    # Polynomial order
        :lsource             => true,
        #---------------------------------------------------------------------------
        #Building matrices
        #---------------------------------------------------------------------------
        :ldss_laplace        => true,
        :lsparse             => true,
        :ldss_differentiation => false,#true, # with true I get error (DSS_generic_matrix not defined)
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :μ                    => [0.1, 0.1],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename         => "./meshes/gmsh_grids/kopriva.msh",
        :gmsh_filename         => "./meshes/gmsh_grids/kopriva_periodic.msh",
        #---------------------------------------------------------------------------
        # grid modification parameters
        #--------------------------------------------------------------------------- 
        :xscale              => 1.0,
        :yscale              => 1.0,
        :xdisp               => 0.0,
        :ydisp               => 0.0,
        #---------------------------------------------------------------------------
        # Refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => true,
        :init_refine_lvl     => 1,
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        #:lwarp               => true,
        #:mount_type          => "agnesi",
        #:a_mount             => 1000.0,
        #:h_mount             => 1.0,
        #:c_mount             => 0.0,
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.15,
        #:mu_y                => 0.15,
        #:filter_type         => "erf",  ##default is erf, use either "erf" for Boyd-Vandeven,"exp" for Warburton Exponential filter, or "quad" for Fischer quadratic filter
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat         => "vtk",
        :loverwrite_output => true,
        :loutput_pert      => true,  #this is only implemented for VTK for now
        :output_dir        => "./output/",
        #:plot_hlines      => [10.0],
        :loutput_pert      => true,
        #---------------------------------------------------------------------------
        # IMEX method
        #---------------------------------------------------------------------------
        :method             => "RK",
        :delta              => 1,
        :k                  => 3,
#        :coeff              => Dict(
#                                   # IMEX Multistep
#                                   :xi       => 2. / 3.,
#                                   :alpha    => alpha,
#                                   :beta     => beta,
#                               ),
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
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
