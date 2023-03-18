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
             neqns, 
             inputs::Dict,
             BCT,
             time,
             T)
    
    dq     = zeros(mesh.npoin, neqns)
    RKcoef = buildRKIntegrator!(TD, T)
    
    for s = 1:length(RKcoef.a)
        
        #
        # rhs[ngl,ngl,nelem]
        #
        rhs_el      = build_rhs(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω,
                                mesh, metrics, T)
        rhs_diff_el = build_rhs_diff(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω,
                                     inputs[:νx], inputs[:νy], mesh, metrics, T)
        #apply_boundary_conditions!(rhs_el, q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)
        #
        # RHS[npoin] = DSS(rhs)
        #
        for ieqn=1:neqns
            RHS = DSSijk_rhs(SD,
                             rhs_el[:,:,:,ieqn] + inputs[:δvisc]*rhs_diff_el[:,:,:,ieqn],
                             mesh.connijk,
                             mesh.nelem, mesh.npoin, mesh.nop,
                             T)
            divive_by_mass_matrix!(RHS, M, QT)
            
            for I=1:mesh.npoin
                dq[I, ieqn] = RKcoef.a[s]*dq[I, ieqn] + Δt*RHS[I]
                q.qn[I, ieqn] = q.qn[I, ieqn] + RKcoef.b[s]*dq[I, ieqn]
            end
            
            #
            # B.C.
            #
        end
        #apply_periodicity!(rhs_el,q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)
        
    end #stages
    
end


function bdf2!(q::St_SolutionVars;
               TD,
               SD,
               QT,
               PT,
               mesh::St_mesh,
               metrics::St_metrics,
               basis, ω,
               M, Δt,
               neqns, 
               inputs::Dict,
               BCT,
               time,
               T)
    
    dq     = zeros(mesh.npoin, neqns)
    RKcoef = buildRKIntegrator!(TD, T)
    rhs_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    
    #
    # rhs[ngl,ngl,nelem]
    #
    #rhs_elnm2 .= rhs_elnm1
    #rhs_elnm2 .= rhs_el
    rhs_eln     = build_rhs(SD, QT, PT, neqns, q.qn,   basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_elnm1   = build_rhs(SD, QT, PT, neqns, q.qnm1, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_elnm2   = build_rhs(SD, QT, PT, neqns, q.qnm2, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    
    for ieq = 1:neqns
        for iel=1:mesh.nelem
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    ip = mesh.connijk[i,j,iel]
                    
                    rhs_el[i,j,iel,ieq] = ω[i]*ω[j]*(18.0*q.qn[ip,ieq] - 9.0*q.qnm1[ip,ieq] + 2.0*q.qnm2[ip,ieq])/11.0 -
                        (6.0*Δt/11.0)*(3.0*rhs_eln[i,j,iel,ieq] - 3.0*rhs_elnm1[i,j,iel,ieq] + rhs_elnm2[i,j,iel,ieq])
                end
            end
        end
    end
                    
    #rhs_diff_el = build_rhs_diff(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω,
    #                             inputs[:νx], inputs[:νy], mesh, metrics, T)
    apply_boundary_conditions!(rhs_el, q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)
    #
    # RHS[npoin] = DSS(rhs)
    #
    for ieqn=1:neqns
        RHS = DSSijk_rhs(SD,
                         rhs_el[:,:,:,ieqn] ,# + inputs[:δvisc]*rhs_diff_el[:,:,:,ieqn],
                         mesh.connijk,
                         mesh.nelem, mesh.npoin, mesh.nop,
                         T)
        divive_by_mass_matrix!(RHS, M, QT)
        
        for I=1:mesh.npoin
            q.qn[I, ieqn] = RHS[I]
        end
        
        #
        # B.C.
        #
    end
    apply_periodicity!(rhs_el,q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)
    
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
                    neqns, 
                    inputs::Dict,
                    BCT,
                    OUTPUT_DIR::String,
                    T)
    it = 0
    t  = inputs[:tinit]
    t0 = t

    plot_at_times = [0.25, 0.5, 1.0, 1.5]    
   
    it_interval = inputs[:diagnostics_interval]
    it_diagnostics = 1
    #
    # RK for first 2 steps
    #
    qp.qnm1 .= qp.qn
    for it = 1:2
        if (mod(it, it_interval) == 0 || it == Nt)
            @printf "   Solution at t = %.6f sec\n" t
            @printf "      min/max(q1) = %.6f %.6f\n" minimum(qp.qn[:,1]) maximum(qp.qn[:,1])
            @printf "      min/max(q2) = %.6f %.6f\n" minimum(qp.qn[:,2]) maximum(qp.qn[:,2])
            @printf "      min/max(q3) = %.6f %.6f\n" minimum(qp.qn[:,3]) maximum(qp.qn[:,3])

            #------------------------------------------
            # Plot initial condition:
            # Notice that I scatter the points to
            # avoid sorting the x and q which would be
            # becessary for a smooth curve plot.
            #------------------------------------------
            title = string( "Tracer: final solution at t=", t)
            jcontour(mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it.", it_diagnostics, ".png"))
            it_diagnostics = it_diagnostics + 1
        end
        t = t0 + Δt
        t0 = t
        
        rk!(qp; TD, SD, QT, PT,
            mesh, metrics, basis, ω, M, Δt, neqns, inputs, BCT, time=t, T)
        
    end
    @printf "----------\n"
    @printf "  n   min/max(qn)  = %.6f %.6f\n" minimum(qp.qn[:,1]) maximum(qp.qn[:,1])
    @printf "  Mn1 min/max(qm1) = %.6f %.6f\n" minimum(qp.qnm1[:,1]) maximum(qp.qnm1[:,1])
    @printf "----------\n"
    #
    # DBF2 from now on:
    #    
    for itbdf = 3:Nt

        #OK don't touch
        qp.qnm2 .= qp.qnm1
        qp.qnm1 .= qp.qn
        #this needs to be here BEFORE the time stepper
        #ok don't touch above.

        rk!(qp; TD, SD, QT, PT,  mesh, metrics, basis, ω, M, Δt, neqns, inputs, BCT, time=t, T)
        #bdf2!(qp; TD, SD, QT, PT, mesh, metrics, basis, ω, M, Δt, neqns, inputs, BCT, time=t, T)
                
        if (mod(itbdf, it_interval) == 0 || itbdf == Nt)
            @printf "   Solution at t = %.6f sec\n" t
            for ieq = 1:neqns
                @printf "      min/max(q[%d]) = %.6f %.6f\n" ieq minimum(qp.qn[:,ieq]) maximum(qp.qn[:,ieq])
            end
            
            #------------------------------------------
            # Plot initial condition:
            # Notice that I scatter the points to
            # avoid sorting the x and q which would be
            # becessary for a smooth curve plot.
            #------------------------------------------
            title = string( "Tracer: final solution at t=", t)
            jcontour(mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it.", it_diagnostics, "N.png"))
            
            it_diagnostics = it_diagnostics + 1
        end
        t = t0 + Δt
        t0 = t
      
    end
      
    #Plot final solution
    title = string( "Tracer: final solution at t=f", inputs[:tend])
    jcontour(mesh.x, mesh.y, qp.qn[:,3], title, string(OUTPUT_DIR, "/END.png"))
    
end
