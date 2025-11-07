function user_inputs()
    # Coefficients of the method
    alpha = Array{Float64, 1}(undef, 2)
    alpha[2] = 4. / 3.
    alpha[1] = - 1. / 3.
    beta = Array{Float64, 1}(undef, 2)
    beta[2] = 2.
    beta[1] = -1.

    alpha = Array{Float64, 1}(undef, 1)
    alpha[1] = 1.
    beta = Array{Float64, 1}(undef, 1)
    beta[1] = 1.

    # Polynomial order
    nop = 4

    # Solver parameters
    solver_par = Dict(:restart  => true,
                      :memory   => 10,
                      :verbose  => 1,
                      :abstol   => 1.e-06,
                      :reltol   => 1.e-06,
                      :maxiters => 100,
                      :prec     => SmoothedAggregationPreconBuilder()
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
    function bcs_fun!(u, time, params)
        SD      = params.SD
        AD      = params.AD
        neqs    = params.neqs
        xmin    = params.mesh.xmin
        xmax    = params.mesh.xmax
        ymin    = params.mesh.ymin
        ymax    = params.mesh.ymax
        zmin    = params.mesh.zmin
        zmax    = params.mesh.zmax

        resetbdyfluxToZero!(params)
        apply_boundary_conditions!(u, params.uaux, time, params.qp.qe,
                                   params.mesh.x, params.mesh.y, params.mesh.z,
                                   params.metrics.nx, params.metrics.ny, params.metrics.nz,
                                   params.mesh.npoin, params.mesh.npoin_linear, 
                                   params.mesh.poin_in_bdy_edge, params.mesh.poin_in_bdy_face,
                                   params.mesh.nedges_bdy, params.mesh.nfaces_bdy, params.mesh.ngl, 
                                   params.mesh.ngr, params.mesh.nelem_semi_inf, params.basis.ψ,
                                   params.basis.dψ, xmax, ymax, zmax, xmin, ymin, zmin, params.RHS,
                                   params.rhs_el, params.ubdy, params.mesh.connijk_lag,
                                   params.mesh.bdy_edge_in_elem, params.mesh.bdy_edge_type,
                                   params.mesh.bdy_face_in_elem, params.mesh.bdy_face_type,
                                   params.mesh.connijk, params.metrics.Jef, params.S_face,
                                   params.S_flux, params.F_surf, params.M_surf_inv, params.M_edge_inv,
                                   params.Minv, params.mp.Tabs, params.mp.qn,
                                   params.ω, neqs, params.inputs, AD, SD)
    end

    # Building fast waves operator
    function build_L!(u, time, params)
        SD      = params.SD
        basis   = params.basis
        ω       = params.ω
        mesh    = params.sem.mesh
        metrics = params.sem.metrics
        N       = nop
        Q       = N + 1
        backend = CPU()

        Le = KernelAbstractions.zeros(backend, TFloat, 1, 1)
        L  = KernelAbstractions.zeros(backend, TFloat, 1, 1)

        Le = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh, metrics, N, Q, TFloat)
        L = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))

        DSS_laplace!(L, SD, Le, ω, mesh, metrics, N, TFloat; llump = true)
    end

    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 4.0, #2π,
        :Δt                   => 0.005,#8.75e-4,
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
        :loverwrite_output => false,
        :loutput_pert      => true,  #this is only implemented for VTK for now
        :output_dir        => "./output/",
        #:plot_hlines      => [10.0],
        :loutput_pert      => true,
        #---------------------------------------------------------------------------
        # IMEX method
        #---------------------------------------------------------------------------
        :method             => "multistep",
        :delta              => 0,
        :k                  => 1,    # add check on k>=1
        :coeff              => Dict(
                                   :xi       => 2. / 3.,
                                   :alpha    => alpha,
                                   :beta     => beta,
                               ),
        :lsolver            => "GMRES",#LinearSolve.KrylovJL_GMRES(),
        :solver_par         => solver_par,
        :S_fun              => S_fun!,
        :L_fun              => L_fun!,
        :bcs_fun            => bcs_fun!,
        :upd_L              => false,
        :build_L            => build_L!,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
