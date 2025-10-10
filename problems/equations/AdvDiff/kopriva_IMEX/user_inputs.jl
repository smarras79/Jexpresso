function user_inputs()
    alpha = Array{Float64, 1}(undef, 2)
    alpha[2] = 4. / 3.
    alpha[1] = - 1. / 3.
    beta = Array{Float64, 1}(undef, 2)
    beta[2] = 2.
    beta[1] = -1.

#    prec = HYPRE.BoomerAMG
    solver_par = (; restart = true,
                    memory = 10,
                    verbose = 1,
                    abstol = 1.e-06,
                    reltol = 1.e-06,
                    maxiters = 100)
#                    N = prec)

    function S_fun!(s_j, u, params, time)
        rhs!(s_j, u, params, time)

        s_j .= params.RHS
    end

    function L_fun!(l_j, u, params, time)
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

    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 1.0, #2π,
        :Δt                   => 0.05,#8.75e-4,
        :ode_solver           => SSPRK54(),
        :diagnostics_at_times => [0.5, 1, 2, 4],
        :output_dir          => "./output/",
        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :nop                 => 4,      # Polynomial order
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
        :loverwrite_output => true,
        :loutput_pert      => true,  #this is only implemented for VTK for now
        :output_dir        => "./output/",
        #:plot_hlines      => [10.0],
        :loutput_pert      => true,
        #---------------------------------------------------------------------------
        # IMEX method
        #---------------------------------------------------------------------------
        :method             => "multistep",
        :delta              => 0,
        :k                  => 2,
        :coeff              => Dict(
                                   :xi       => 2. / 3.,
                                   :alpha    => alpha,
                                   :beta     => beta,
                               ),
        :lsolver            => "GMRES",#LinearSolve.KrylovJL_GMRES(),
        :solver_par         => solver_par,
        :S_fun              => S_fun,
        :L_fun              => L_fun,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
