include("../abstractTypes.jl")
include("../infrastructure/element_matrices.jl")

include("../../io/plotting/jeplots.jl")

mutable struct RK_Integrator{TFloat}
  a::Array{TFloat}
  b::Array{TFloat}
  c::Array{TFloat}
end

function buildRKIntegrator(RKtype::RK3, TFloat)

    RKcoef = RK_Integrator{TFloat}(zeros(TFloat,3),zeros(TFloat,3),zeros(TFloat,3))
    
    RKcoef.a=[0.0 -5/9 -153/128]
    RKcoef.b=[0.0 1/3 3/4]
    RKcoef.c=[1/3 15/16 8/15]
    
    return RKcoef
end


function buildRKIntegrator!(RKtype::RK5, TFloat)

    RKcoef = RK_Integrator{TFloat}(zeros(TFloat,5),zeros(TFloat,5),zeros(TFloat,5))
    
    RKcoef.a = [TFloat(0),
                TFloat(-567301805773)  / TFloat(1357537059087),
                TFloat(-2404267990393) / TFloat(2016746695238),
                TFloat(-3550918686646) / TFloat(2091501179385),
                TFloat(-1275806237668) / TFloat(842570457699 )]

    RKcoef.b = [TFloat(1432997174477) / TFloat(9575080441755 ),
                TFloat(5161836677717) / TFloat(13612068292357),
                TFloat(1720146321549) / TFloat(2090206949498 ),
                TFloat(3134564353537) / TFloat(4481467310338 ),
                TFloat(2277821191437) / TFloat(14882151754819)]

    RKcoef.c = [TFloat(0),
                TFloat(1432997174477) / TFloat(9575080441755),
                TFloat(2526269341429) / TFloat(6820363962896),
                TFloat(2006345519317) / TFloat(3224310063776),
                TFloat(2802321613138) / TFloat(2924317926251)]

    return RKcoef

end

function rk!(q::St_SolutionVars;
             TD,
             SD,
             QT,
             PT,
             mesh::St_mesh,
             metrics::St_metrics,
             basis, ω,
             M, Δt,
             nvars, 
             inputs::Dict,
             T)
    
    dq     = zeros(mesh.npoin, nvars)
    RKcoef = buildRKIntegrator!(TD, T)
    
    for s = 1:length(RKcoef.a)
        
        #
        # rhs[ngl,ngl,nelem]
        #
        rhs_el      = build_rhs(SD, QT, PT, nvars, q, basis.ψ, basis.dψ, ω, mesh, metrics, T)
        rhs_diff_el = build_rhs_diff(SD, QT, PT, nvars, q, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
        
        #
        # RHS[npoin] = DSS(rhs)
        #
        for ivar=1:nvars
            RHS = DSSijk_rhs(SD,
                             rhs_el[:,:,:,ivar] + rhs_diff_el[:,:,:,ivar],
                             mesh.connijk,
                             mesh.nelem, mesh.npoin, mesh.nop,
                             T)
            divive_by_mass_matrix!(RHS, M, QT)
            
            for I=1:mesh.npoin
                dq[I, ivar] = RKcoef.a[s]*dq[I, ivar] + Δt*RHS[I]
                q.qn[I, ivar] = q.qn[I, ivar] + RKcoef.b[s]*dq[I, ivar]
            end
            
            #
            # B.C.
            #
        end
        apply_boundary_conditions!(q, mesh, inputs, SD)
        
        end #stages

    #return qp
    
end


function rk4!(q::St_SolutionVars;
             TD,
             SD,
             QT,
             PT,
             mesh::St_mesh,
             metrics::St_metrics,
             basis, ω,
             M, Δt,
             nvars, 
             inputs::Dict,
             T)
    
    dq      = zeros(T, mesh.npoin, nvars)
    Q1 = Q2 = zeros(T, mesh.npoin, nvars)

    #
    # Stage 1
    #   
    rhs_el      = build_rhs(SD, QT, PT, nvars, q, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, nvars, q, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
    for ivar=1:nvars
        RHS = DSSijk_rhs(SD,
                         rhs_el[:,:,:,ivar] + rhs_diff_el[:,:,:,ivar],
                         mesh.connijk,
                         mesh.nelem, mesh.npoin, mesh.nop,
                         T)
        divive_by_mass_matrix!(RHS, M, QT)
        
        Q1[I,ivar] = q.qn[I,ivar] + Δt*RHS[I]
        apply_boundary_conditions!(Q1[:,ivar], mesh, inputs, SD)
    end

    #
    # Stage 2
    #   
    rhs_el      = build_rhs(SD, QT, PT, nvars, Q1, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, nvars, Q1, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
    for ivar=1:nvars
        RHS = DSSijk_rhs(SD,
                         rhs_el[:,:,:,ivar] + rhs_diff_el[:,:,:,ivar],
                         mesh.connijk,
                         mesh.nelem, mesh.npoin, mesh.nop,
                         T)
        divive_by_mass_matrix!(RHS, M, QT)
        
        Q2[:,ivar] = 0.75*q.qn[:,ivar] + 0.25*(Q1[:, ivar] + Δt*RHS[I])
        apply_boundary_conditions!(Q2[:,ivar], mesh, inputs, SD)
    end

    
    #
    # Stage 3: qⁿ⁺¹
    #   
    rhs_el      = build_rhs(SD, QT, PT, nvars, Q1, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, nvars, Q1, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
    for ivar=1:nvars
        RHS = DSSijk_rhs(SD,
                         rhs_el[:,:,:,ivar] + rhs_diff_el[:,:,:,ivar],
                         mesh.connijk,
                         mesh.nelem, mesh.npoin, mesh.nop,
                         T)
        divive_by_mass_matrix!(RHS, M, QT)
        
        q.qn[:,ivar] = (1.0/3.0)*q.qn[:,ivar] + (2.0/3.0)*(Q2[:, ivar] + Δt*RHS[I])
        apply_boundary_conditions!(q.qn[:,ivar], mesh, inputs, SD)
    end
    
    
end

function time_loop!(TD,
                    SD,
                    QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp,
                    M,
                    Nt, Δt,
                    nvars, 
                    inputs::Dict,
                    OUTPUT_DIR::String,
                     T)
    it = 0
    t  = inputs[:tinit]
    t0 = t

    plot_at_times = [0.25, 0.5, 1.0, 1.5]    
   
    it_interval = inputs[:diagnostics_interval]
    it_diagnostics = 1
    for it = 1:Nt 

        qp.qnm1[:,:] = qp.qn[:,:] #[npoin, nvar]
        rk!(qp; TD, SD, QT, PT, mesh, metrics, basis, ω, M, Δt, nvars, inputs, T)
        #rk4!(qp; TD, SD, QT, PT, mesh, metrics, basis, ω, M, Δt, nvars, inputs, T)
       
        if (mod(it, it_interval) == 0 || it == Nt)
            @printf "   Solution at t = %.6f sec\n" t
            @printf "      min(qⁿ), max(qⁿ) = %.6f, %.6f\n" minimum(qp.qn[:,1]) maximum(qp.qn[:,1])
            
            #------------------------------------------
            # Plot initial condition:
            # Notice that I scatter the points to
            # avoid sorting the x and q which would be
            # becessary for a smooth curve plot.
            #------------------------------------------
            title = string( "Tracer: final solution at t=", t)
            jcontour(mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it.", it_diagnostics, ".png"))
            #title = string( "Tracer: final solution at t=", t)
            #jcontour(mesh.x, mesh.y, qp.qnm1[:,1], title, string(OUTPUT_DIR, "/it.", it_diagnostics, "Nm1.png"))

            it_diagnostics = it_diagnostics + 1
        end
        t = t0 + Δt
        t0 = t
       
    end
      
    #Plot final solution
    title = string( "Tracer: final solution at t=%.8f", inputs[:tend])
    jcontour(mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/END.png"))
    
end
