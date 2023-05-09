#---------------------------------------------------------------------------
# Fetch equations name to access the user_rhs functions
#---------------------------------------------------------------------------
if (length(ARGS) === 1) #equations
    user_flux_dir   = string("../../equations/", ARGS[1], "/user_flux.jl")
    user_source_dir = string("../../equations/", ARGS[1], "/user_source.jl")
elseif (length(ARGS) === 2)  #equations/equations_case_name
    user_flux_dir   = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_flux.jl")
    user_source_dir = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl")
end
include(user_flux_dir)
include(user_source_dir)
include("../ArtificialViscosity/DynSGS.jl")
#---------------------------------------------------------------------------

function rhs!(du, u, params, time)
    
    RHS = build_rhs(params.SD, params.QT, params.PT,
                    u,
                    params.neqs,
                    params.basis, params.ω,
                    params.mesh, params.metrics,
                    params.M, params.De, params.Le,
                    time,
                    params.inputs, params.Δt, params.deps, params.T;
                    qnm1=params.qnm1, qnm2=params.qnm2)
    for i=1:params.neqs
       idx = (i-1)*params.mesh.npoin
       du[idx+1:i*params.mesh.npoin] = @view RHS[:,i]
    end  
    return du #This is already DSSed
end

#--------------------------------------------------------------------------------------------------------------------------------------------------
# AdvDiff
#--------------------------------------------------------------------------------------------------------------------------------------------------
function build_rhs_matrix_formulation(SD::NSD_1D, QT::Inexact, PT::AdvDiff, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T)
    
    Fuser = user_flux(T, SD, qp, mesh)
    
    #
    # Linear RHS in flux form: f = u*u
    #  
    RHS = zeros(mesh.npoin)
    qe  = zeros(mesh.ngl)
    fe  = zeros(mesh.ngl)
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            qe[i] = qp[I,1]
            fe[i] = Fuser[I]
        end
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            for j=1:mesh.ngl
                RHS[I] = RHS[I] - De[i,j,iel]*fe[j] + inputs[:νx]*Le[i,j,iel]*qe[j]
            end
        end
    end
    
    # M⁻¹*rhs where M is diagonal
    RHS .= RHS./M

    apply_periodicity!(SD, RHS, qp, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, 0, neqs)
    
    return RHS
    
end

function build_rhs(SD::NSD_1D, QT::Inexact, PT::AdvDiff, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T)

    F      = zeros(mesh.ngl, mesh.nelem, neqs)
    dFdξ   = zeros(neqs)
    rhs_el = zeros(mesh.ngl, mesh.nelem, neqs)
    qq     = zeros(mesh.npoin, neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= qp[idx+1:i*mesh.npoin]
    end    
    Fuser = user_flux(T, SD, qq, mesh; neqs=neqs)

    for iel=1:mesh.nelem
        
        Jac     = mesh.Δx[iel]/2
        metrics = 2.0/mesh.Δx[iel]
        
        for i=1:mesh.ngl
            ip = mesh.conn[i,iel]
            F[i,iel,1:neqs] .= Fuser[ip,1:neqs]
        end

        for i=1:mesh.ngl
            dFdξ[1:neqs] .= 0.0
            for k = 1:mesh.ngl
                dFdξ[1:neqs] .= dFdξ[1:neqs] .+ basis.dψ[k,i]*F[k, iel, 1:neqs]*metrics
            end
            rhs_el[i, iel, 1:neqs] .+= ω[i]*Jac*dFdξ[1:neqs]
        end
    end
    
    RHS = DSS_rhs(SD, rhs_el, mesh.conn, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    divive_by_mass_matrix!(RHS, M, QT, neqs)
    
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qp[idx+1:i*mesh.npoin] .= qq[:,i]
    end
    apply_periodicity!(SD, RHS, qp, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, 0, neqs)
    
    return RHS
end

function build_rhs(SD::NSD_1D, QT::Exact, PT::AdvDiff, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T) nothing end

function build_rhs(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T)
    
    F      = zeros(T, mesh.ngl, mesh.ngl, mesh.nelem)
    G      = zeros(T, mesh.ngl, mesh.ngl, mesh.nelem)
    rhs_el = zeros(T, mesh.ngl, mesh.ngl, mesh.nelem)
    
    #B.C.
    apply_boundary_conditions!(SD, rhs_el, qp, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, time, neqs)   
    Fuser, Guser = user_flux(T, SD, qp, mesh)
    
    for iel=1:mesh.nelem
        for j=1:mesh.ngl, i=1:mesh.ngl
            ip = mesh.connijk[i,j,iel]
            
            F[i,j,iel] = Fuser[ip]
            G[i,j,iel] = Guser[ip]
        end
    end
    
    # for ieq = 1:neqs
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                for k = 1:mesh.ngl
                    dFdξ = dFdξ + basis.dψ[k, i]*F[k,j,iel]
                    dFdη = dFdη + basis.dψ[k, j]*F[i,k,iel]

                    dGdξ = dGdξ + basis.dψ[k, i]*G[k,j,iel]
                    dGdη = dGdη + basis.dψ[k, j]*G[i,k,iel]
                end
                dFdx = dFdξ*metrics.dξdx[i,j,iel] + dFdη*metrics.dηdx[i,j,iel]
                dGdy = dGdξ*metrics.dξdy[i,j,iel] + dGdη*metrics.dηdy[i,j,iel]
                rhs_el[i, j, iel] -= ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx + dGdy)
            end
        end
    end
    #end
    #show(stdout, "text/plain", el_matrices.D)
    RHS = DSS_rhs(SD, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin,neqs, mesh.nop, T)
    
    #Build rhs_el(diffusion)
    if (inputs[:lvisc] == true)
        μ = zeros(mesh.nelem,1)
        if (inputs[:visc_model] === "dsgs")
            if (rem(time, Δt) < 5e-4 && time > 0.0)
                global  q2 .= q1
                global  q1 .= q3
                global  q3 .= qq
            end
            compute_viscosity!(mu, SD, PT, q3, q1, q2, RHS, Δt, mesh, metrics)
        else
            μ[:] .= inputs[:νx]
        end
        rhs_diff_el = build_rhs_diff(SD, QT, PT, qp,  neqs, basis, ω, inputs, mesh, metrics, μ, T;)
        RHS .= RHS .+ DSS_rhs(SD, rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    end
    
    divive_by_mass_matrix!(RHS, M, QT,neqs)
    
    return RHS
    
end

function build_rhs(SD::NSD_2D, QT::Exact, PT::AdvDiff, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T) nothing end

#--------------------------------------------------------------------------------------------------------------------------------------------------
# LinearCLaw
#--------------------------------------------------------------------------------------------------------------------------------------------------
function build_rhs(SD::NSD_2D, QT::Inexact, PT::LinearCLaw, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T)    
    
    F    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    G    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)

    rhs_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    qq = zeros(mesh.npoin,neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= qp[idx+1:i*mesh.npoin]
    end
    Fuser, Guser = user_flux(T, SD, qq, mesh)
    dFdx = zeros(neqs)
    dFdξ = zeros(neqs)
    dGdξ = zeros(neqs)
    dGdy = zeros(neqs)
    dFdη = zeros(neqs)
    dGdη = zeros(neqs)
    
    for iel=1:mesh.nelem

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]
                F[i,j,iel,1] = Fuser[ip,1]
                F[i,j,iel,2] = Fuser[ip,2]
                F[i,j,iel,3] = Fuser[ip,3]

                G[i,j,iel,1] = Guser[ip,1]
                G[i,j,iel,2] = Guser[ip,2]
                G[i,j,iel,3] = Guser[ip,3]

            end
        end

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                dFdξ = zeros(T, neqs)
                dFdη = zeros(T, neqs)
                dGdξ = zeros(T, neqs) 
                dGdη = zeros(T, neqs)
                for k = 1:mesh.ngl
                    dFdξ[1:neqs] .= dFdξ[1:neqs] .+ basis.dψ[k,i]*F[k,j,iel,1:neqs]
                    dFdη[1:neqs] .= dFdη[1:neqs] .+ basis.dψ[k,j]*F[i,k,iel,1:neqs]
                    
                    dGdξ[1:neqs] .= dGdξ[1:neqs] .+ basis.dψ[k,i]*G[k,j,iel,1:neqs]
                    dGdη[1:neqs] .= dGdη[1:neqs] .+ basis.dψ[k,j]*G[i,k,iel,1:neqs]
                end
                dFdx .= dFdξ[1:neqs]*metrics.dξdx[i,j,iel] .+ dFdη[1:neqs]*metrics.dηdx[i,j,iel]
                dGdy .= dGdξ[1:neqs]*metrics.dξdy[i,j,iel] .+ dGdη[1:neqs]*metrics.dηdy[i,j,iel]
                rhs_el[i,j,iel,1:neqs] .-= ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx[1:neqs] .+ dGdy[1:neqs])
            end
        end
    end

    apply_boundary_conditions!(SD, rhs_el, qq, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, Δt*(floor(time/Δt)), neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qp[idx+1:i*mesh.npoin] .= qq[:,i]
    end
    RHS = DSS_rhs(SD, rhs_el .+ rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    
    if (inputs[:lvisc] == true)
        μ = zeros(mesh.nelem,1)
        if (inputs[:visc_model] === "dsgs")
            compute_viscosity!(mu, SD, PT, q3, q1, q2, RHS, Δt, mesh, metrics)
        else
            μ[:] .= inputs[:νx]
        end
        if (rem(time, Δt) < 5e-4 && time > 0.0)
            global  q2 .= q1
            global  q1 .= q3
            global  q3 .= qq
        end
        rhs_diff_el = build_rhs_diff(SD, QT, PT, qp,  neqs, basis, ω, inputs, mesh, metrics, μ, T;)
        RHS .= RHS .+ DSS_rhs(SD, rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    end
    
    divive_by_mass_matrix!(RHS, M, QT,neqs)
    
    return RHS
end


#--------------------------------------------------------------------------------------------------------------------------------------------------
# ShallowWater
#--------------------------------------------------------------------------------------------------------------------------------------------------
function build_rhs(SD::NSD_1D, QT::Inexact, PT::ShallowWater, qp::Array, neqs, basis, ω,
                   mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T)

    F      = zeros(mesh.ngl,mesh.nelem, neqs)
    F1     = zeros(mesh.ngl,mesh.nelem, neqs)
    rhs_el = zeros(mesh.ngl,mesh.nelem, neqs)
    qq     = zeros(mesh.npoin,neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= qp[idx+1:i*mesh.npoin]
    end
    qq[:,1] = max.(qq[:,1],0.001)
    qq[:,2] = max.(qq[:,2],0.0)
    #S =  user_source_friction(SD, T, qq, mesh.npoin)
   
    Fuser, Fuser1 = user_flux(T, SD, qq, mesh)
    dFdx = zeros(neqs)
    dFdξ = zeros(neqs)
    gHsx = zeros(neqs)
    for iel=1:mesh.nelem
        dξdx = 2.0/mesh.Δx[iel]
        for i=1:mesh.ngl
                ip = mesh.conn[i,iel]
                F[i,iel,1] = Fuser[ip,1]
                F[i,iel,2] = Fuser[ip,2]

                F1[i,iel,1] = Fuser1[ip,1]
                F1[i,iel,2] = Fuser1[ip,2]
                #@info Fuser[ip,1] + Fuser1[ip,1], Fuser[ip,2] + Fuser1[ip,2]
        end
        for i=1:mesh.ngl
            dFdξ = zeros(T, neqs)
            dFdξ1 = zeros(T, neqs)
            for k = 1:mesh.ngl
                dFdξ[1:neqs] .= dFdξ[1:neqs] .+ basis.dψ[k,i]*F[k,iel,1:neqs]*dξdx

                dFdξ1[1:neqs] .= dFdξ1[1:neqs] .+ basis.dψ[k,i]*F1[k,iel,1:neqs]*dξdx
                #@info i,dFdξ[1:neqs], dFdξ1[1:neqs]
            end
            ip = mesh.conn[i,iel]
            x = mesh.x[ip]
            Hb = zb[ip]
            Hs = max(qq[ip,1] - Hb,0.001)
            gHsx[1] = 1.0
            gHsx[2] = qq[ip,1]*9.81#Hs*9.81
            dFdx[1] = gHsx[1] * (dFdξ[1]) + dFdξ1[1] 
            dFdx[2] = gHsx[2] * (dFdξ[2]) + dFdξ1[2] #+ S[ip]*qq[ip,1]*9.81
            rhs_el[i,iel,1:neqs] .-= ω[i]*mesh.Δx[iel]/2*dFdx[1:neqs]
        end
    end   
    apply_boundary_conditions!(SD, rhs_el, qq, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, Δt*(floor(time/Δt)), neqs)
    RHS = DSS_rhs(SD, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)

    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qp[idx+1:i*mesh.npoin] .= qq[:,i]
    end
    
    if (inputs[:lvisc] == true)
        mu = zeros(mesh.nelem,1)
        if (inputs[:visc_model] === "dsgs")
            compute_viscosity!(mu, SD, PT, q3, q1, q2, RHS, Δt, mesh, metrics)
        else
            mu[:] = inputs[:νx]
        end
        if (rem(time, Δt) < 5e-4 && time > 0.0)
            global  q2 .= q1
            global  q1 .= q3
            global  q3 .= qq
        end
        rhs_diff_el = build_rhs_diff(SD, QT, PT, qp,  neqs, basis, ω, inputs, mesh, metrics, mu, T;)
        RHS .= RHS .+ DSS_rhs(SD, rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    end
    
    divive_by_mass_matrix!(RHS, M, QT,neqs)

    return RHS
end

function build_rhs(SD::NSD_2D, QT::Inexact, PT::ShallowWater, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T)
    F    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    G    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    F1    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    G1    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    rhs_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    qq = zeros(mesh.npoin,neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= qp[idx+1:i*mesh.npoin]
    end
    Fuser, Guser, Fuser1, Guser1 = user_flux(T, SD, qq, mesh)
    dFdx = zeros(neqs)
    dFdξ = zeros(neqs)
    dGdξ = zeros(neqs)
    dGdy = zeros(neqs)
    dFdη = zeros(neqs)
    dGdη = zeros(neqs)
    gHsx = zeros(neqs)
    gHsy = zeros(neqs)
    for iel=1:mesh.nelem

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]
                F[i,j,iel,1] = Fuser[ip,1]
                F[i,j,iel,2] = Fuser[ip,2]
                F[i,j,iel,3] = Fuser[ip,3]
                
                F1[i,j,iel,1] = Fuser1[ip,1] 
                F1[i,j,iel,2] = Fuser1[ip,2] 
                F1[i,j,iel,3] = Fuser1[ip,3]  
               
                G[i,j,iel,1] = Guser[ip,1]
                G[i,j,iel,2] = Guser[ip,2]
                G[i,j,iel,3] = Guser[ip,3]

                G1[i,j,iel,1] = Guser1[ip,1] 
                G1[i,j,iel,2] = Guser1[ip,2] 
                G1[i,j,iel,3] = Guser1[ip,3] 
            end
        end

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                dFdξ = zeros(T, neqs)
                dFdη = zeros(T, neqs)
                dGdξ = zeros(T, neqs)
                dGdη = zeros(T, neqs)
                dFdξ1 = zeros(T, neqs)
                dFdη1 = zeros(T, neqs)
                dGdξ1 = zeros(T, neqs)
                dGdη1 = zeros(T, neqs)
                for k = 1:mesh.ngl
                    dFdξ[1:neqs] .= dFdξ[1:neqs] .+ basis.dψ[k,i]*F[k,j,iel,1:neqs]
                    dFdη[1:neqs] .= dFdη[1:neqs] .+ basis.dψ[k,j]*F[i,k,iel,1:neqs]
                     
                    dFdξ1[1:neqs] .= dFdξ1[1:neqs] .+ basis.dψ[k,i]*F1[k,j,iel,1:neqs]
                    dFdη1[1:neqs] .= dFdη1[1:neqs] .+ basis.dψ[k,j]*F1[i,k,iel,1:neqs]    

                    dGdξ[1:neqs] .= dGdξ[1:neqs] .+ basis.dψ[k,i]*G[k,j,iel,1:neqs]
                    dGdη[1:neqs] .= dGdη[1:neqs] .+ basis.dψ[k,j]*G[i,k,iel,1:neqs]

                    dGdξ1[1:neqs] .= dGdξ1[1:neqs] .+ basis.dψ[k,i]*G1[k,j,iel,1:neqs]
                    dGdη1[1:neqs] .= dGdη1[1:neqs] .+ basis.dψ[k,j]*G1[i,k,iel,1:neqs]
                end
                ip = mesh.connijk[i,j,iel]
                x = mesh.x[ip]
                y = mesh.y[ip]
                Hb = bathymetry(x,y)
                Hs = max(qq[ip,1] - Hb,0.001)
                gHsx[1] = 1.0
                gHsx[2] = Hs * 9.81
                gHsx[3] = 1.0
                gHsy[1] = 1.0
                gHsy[2] = 1.0
                gHsy[3] = Hs * 9.81
                dFdx .= gHsx .* (dFdξ[1:neqs]*metrics.dξdx[i,j,iel] .+ dFdη[1:neqs]*metrics.dηdx[i,j,iel]) + dFdξ1[1:neqs]*metrics.dξdx[i,j,iel] .+ dFdη1[1:neqs]*metrics.dηdx[i,j,iel]
                dGdy .= gHsy .* (dGdξ[1:neqs]*metrics.dξdy[i,j,iel] .+ dGdη[1:neqs]*metrics.dηdy[i,j,iel]) + dGdξ1[1:neqs]*metrics.dξdy[i,j,iel] .+ dGdη1[1:neqs]*metrics.dηdy[i,j,iel]
                rhs_el[i,j,iel,1:neqs] .-= ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx[1:neqs] .+ dGdy[1:neqs])
            end
        end
    end

    if (inputs[:lvisc] == true)
        μ = zeros(mesh.nelem,1)
        
        if (lowercase(inputs[:visc_model]) === "dsgs")
            
            if (rem(time, Δt) == 0 && time > 0.0)
                global  q1 .= q2
                global  q2 .= qq
            end
            
            compute_viscosity!(μ, SD, PT, qq, q1, q2, RHS, Δt, mesh, metrics)
        else
            μ[:] .= inputs[:νx]
        end

        rhs_diff_el = build_rhs_diff(SD, QT, PT, qp, neqs, basis, ω, inputs, mesh, metrics, μ, T;)
        RHS .= RHS .+ DSS_rhs(SD, rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    end
    

    apply_boundary_conditions!(SD, rhs_el, qq, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, Δt*(floor(time/Δt)), neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qp[idx+1:i*mesh.npoin] .= qq[:,i]
    end
    RHS = DSS_rhs(SD, rhs_el .+ rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    divive_by_mass_matrix!(RHS, M, QT,neqs)
    return RHS
end


function build_rhs(SD::NSD_1D, QT::Inexact, PT::SoilTopo, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T)
    F      = zeros(mesh.ngl,mesh.nelem)
    F1     = zeros(mesh.ngl,mesh.nelem)
    rhs_el = zeros(mesh.ngl,mesh.nelem)
    dFdx = 0.0
    dFdξ = 0.0
    dFdξ1 = 0.0
    S =  user_source_friction(SD, T, deps, mesh.npoin)
    for iel=1:mesh.nelem
        dξdx = 2.0/mesh.Δx[iel]
        for i=1:mesh.ngl
                ip = mesh.conn[i,iel]
                F[i,iel] = deps[ip,1]
                F1[i,iel] = zb[ip]
        end
        for i=1:mesh.ngl
                dFdξ = 0.0
                dFdξ1 = 0.0
                for k = 1:mesh.ngl
                    dFdξ = dFdξ + basis.dψ[k,i]*F[k,iel]*dξdx
                    dFdξ1 = dFdξ1 + basis.dψ[k,i] * F1[k,iel]*dξdx
                end
                ip = mesh.conn[i,iel]
                x = mesh.x[ip]
                #if (deps[ip,1] > 0.05)
                    factor = (deps[ip,2]^2/(9.81*deps[ip,1]^3+1e-16)-1)
                #else
                 #   factor = 0.0
                #end
                dFdx = dFdξ1 - factor * dFdξ + S[ip]
                rhs_el[i,iel] += ω[i]*mesh.Δx[iel]/2*dFdx
        end
    end
    
    RHS = DSS_rhs(SD, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    divive_by_mass_matrix!(RHS, M, QT,neqs)
    
    return RHS
end

#--------------------------------------------------------------------------------------------------------------------------------------------------
# CompEuler:
#--------------------------------------------------------------------------------------------------------------------------------------------------
function build_rhs(SD::NSD_1D, QT::Inexact, PT::CompEuler, qp::Array, neqs, basis, ω,
                   mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T)

    F      = zeros(T, mesh.ngl,mesh.nelem, neqs)
    rhs_el = zeros(T, mesh.ngl,mesh.nelem, neqs)
    qq     = zeros(T, mesh.npoin,neqs)
    
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= qp[idx+1:i*mesh.npoin]
    end
    Fuser = user_flux(T, SD, qq, mesh; neqs=neqs)
    
    dFdξ = zeros(neqs)
    for iel=1:mesh.nelem
        dξdx = 2.0/mesh.Δx[iel]
        for i=1:mesh.ngl
            ip = mesh.conn[i,iel]
            F[i,iel,1] = Fuser[ip,1]
            F[i,iel,2] = Fuser[ip,2]
            F[i,iel,3] = Fuser[ip,3]
        end
        for i=1:mesh.ngl
            dFdξ = zeros(T, neqs)
            for k = 1:mesh.ngl
                dFdξ[1:neqs] .= dFdξ[1:neqs] .+ basis.dψ[k,i]*F[k,iel,1:neqs]*dξdx 
            end
            
            rhs_el[i,iel,1:neqs] .-= ω[i]*mesh.Δx[iel]/2*dFdξ[1:neqs]
        end
    end
    
    apply_boundary_conditions!(SD, rhs_el, qq, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, Δt*(floor(time/Δt)), neqs)
    RHS = DSS_rhs(SD, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)

    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qp[idx+1:i*mesh.npoin] .= qq[:,i]
    end
    if (inputs[:lvisc] == true)
        μ = zeros(mesh.nelem,1)
        
        if (inputs[:visc_model] === "dsgs")
            
            if (rem(time, Δt) == 0 && time > 0.0)
                #global  q1 .= q2
                #global  q2 .= qq
                qnm1 .= qnm2
                qnm2 .= qq
            end
            
            μ = compute_viscosity(SD, PT, qq, q1, q2, RHS, Δt, mesh, metrics)
        else
            μ[:] .= inputs[:νx]
        end
        rhs_diff_el = build_rhs_diff(SD, QT, PT, qp, neqs, basis, ω, inputs, mesh, metrics, μ, T;)
        RHS .= RHS .+ DSS_rhs(SD, rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    end
    divive_by_mass_matrix!(RHS, M, QT,neqs)
    
    return RHS
end


### OLD
function oldbuild_rhs(SD::NSD_2D, QT::Inexact, PT::CompEuler, qp::Array, neqs, basis, ω,
                   mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T; qnm1=zeros(1,1), qnm2=zeros(1,1))

    F      = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    G      = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    S      = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    rhs_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    qq = zeros(mesh.npoin,neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= qp[idx+1:i*mesh.npoin]
    end
    Fuser, Guser = olduser_flux(T, SD, qq, mesh; neqs=neqs)
    Suser        = olduser_source(T, qq, mesh.npoin; neqs=neqs)
    
    dFdx = zeros(neqs)
    dGdy = zeros(neqs)
    dFdξ = zeros(neqs)
    dGdξ = zeros(neqs)
    dFdη = zeros(neqs)
    dGdη = zeros(neqs)
    for iel=1:mesh.nelem

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]
                F[i,j,iel,1:neqs] = Fuser[ip,1:neqs]
                G[i,j,iel,1:neqs] = Guser[ip,1:neqs]
                S[i,j,iel,1:neqs] = Suser[ip,1:neqs]
            end
        end

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                dFdξ = zeros(T, neqs)
                dFdη = zeros(T, neqs)
                dGdξ = zeros(T, neqs)
                dGdη = zeros(T, neqs)
                for k = 1:mesh.ngl
                    dFdξ[1:neqs] .= dFdξ[1:neqs] .+ basis.dψ[k,i]*F[k,j,iel,1:neqs]
                    dFdη[1:neqs] .= dFdη[1:neqs] .+ basis.dψ[k,j]*F[i,k,iel,1:neqs]
                    
                    dGdξ[1:neqs] .= dGdξ[1:neqs] .+ basis.dψ[k,i]*G[k,j,iel,1:neqs]
                    dGdη[1:neqs] .= dGdη[1:neqs] .+ basis.dψ[k,j]*G[i,k,iel,1:neqs]
                end
                
                dFdx .= dFdξ[1:neqs]*metrics.dξdx[i,j,iel] .+ dFdη[1:neqs]*metrics.dηdx[i,j,iel]
                dGdy .= dGdξ[1:neqs]*metrics.dξdy[i,j,iel] .+ dGdη[1:neqs]*metrics.dηdy[i,j,iel]
                rhs_el[i,j,iel,1:neqs] .-= ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx[1:neqs] .+ dGdy[1:neqs]) .-=ω[i]*ω[j]*metrics.Je[i,j,iel]*S[i,j,iel,1:neqs] #gravity
            end
        end
    end
    
    apply_boundary_conditions!(SD, rhs_el, qq, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, Δt*(floor(time/Δt)), neqs)
    RHS = DSS_rhs(SD, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)

    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qp[idx+1:i*mesh.npoin] .= qq[:,i]
    end
    
    if (inputs[:lvisc] == true)
        μ = zeros(mesh.nelem,1)

        if (lowercase(inputs[:visc_model]) === "dsgs")
            
            if (rem(time, Δt) == 0 && time > 0.0)
                global  q1 .= q2
                global  q2 .= qq
            end
            
            compute_viscosity!(μ, SD, PT, qq, q1, q2, RHS, Δt, mesh, metrics)
        else
            μ[:] .= inputs[:νx]
        end
        rhs_diff_el = build_rhs_diff(SD, QT, PT, qp, neqs, basis, ω, inputs, mesh, metrics, μ, T;)
        RHS .= RHS .+ DSS_rhs(SD, rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    end
    
    divive_by_mass_matrix!(RHS, M, QT,neqs)
    
    return RHS
end

### OPTIMAL
function build_rhs(SD::NSD_2D, QT::Inexact, PT::CompEuler, qp::Array, neqs, basis, ω,
                   mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T; qnm1=zeros(1,1), qnm2=zeros(1,1))

    F      = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    G      = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    S      = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    rhs_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    qq = zeros(mesh.npoin,neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= qp[idx+1:i*mesh.npoin]
    end
    
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            ip = mesh.connijk[i,j,iel]
            
            F[i,j,iel,1:neqs], G[i,j,iel,1:neqs] = user_flux(T, SD, qq[ip,1:neqs], mesh; neqs=neqs)
            S[i,j,iel,1:neqs] = user_source(T, qq[ip,1:neqs], mesh.npoin; neqs=neqs)
        end

        for ieq = 1:neqs
            for j=1:mesh.ngl, i=1:mesh.ngl
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                for k = 1:mesh.ngl
                    dFdξ += basis.dψ[k,i]*F[k,j,iel,ieq]
                    dFdη += basis.dψ[k,j]*F[i,k,iel,ieq]
                    
                    dGdξ += basis.dψ[k,i]*G[k,j,iel,ieq]
                    dGdη += basis.dψ[k,j]*G[i,k,iel,ieq]
                end

                dFdx = dFdξ*metrics.dξdx[i,j,iel] + dFdη*metrics.dηdx[i,j,iel]
                dGdy = dGdξ*metrics.dξdy[i,j,iel] + dGdη*metrics.dηdy[i,j,iel]
                rhs_el[i,j,iel,ieq] -= ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx + dGdy) -+ ω[i]*ω[j]*metrics.Je[i,j,iel]*S[i,j,iel,ieq] #gravity
            end
        end
    end
    
    apply_boundary_conditions!(SD, rhs_el, qq, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, Δt*(floor(time/Δt)), neqs)
    RHS = DSS_rhs(SD, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)

    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qp[idx+1:i*mesh.npoin] .= qq[:,i]
    end
    
    if (inputs[:lvisc] == true)
        μ = zeros(mesh.nelem,1)

        if (lowercase(inputs[:visc_model]) === "dsgs")
            
            if (rem(time, Δt) == 0 && time > 0.0)
                qnm1 .= qnm2
                qnm2 .= qq
            end
            
            compute_viscosity!(μ, SD, PT, qq, qnm1, qnm2, RHS, Δt, mesh, metrics)
        else
            μ[:] .= inputs[:νx]
        end
        rhs_diff_el = build_rhs_diff(SD, QT, PT, qp, neqs, basis, ω, inputs, mesh, metrics, μ, T;)
        RHS .= RHS .+ DSS_rhs(SD, rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    end
    
    divive_by_mass_matrix!(RHS, M, QT,neqs)
    
    return RHS
end


#--------------------------------------------------------------------------------------------------------------------------------------------------
# Source terms:
#--------------------------------------------------------------------------------------------------------------------------------------------------
function build_rhs_source(SD::NSD_2D,
                          QT::Inexact,
                          q::Array,
                          mesh::St_mesh,
                          M::AbstractArray, #M is sparse for exact integration
                          T)

    S = user_source(q, mesh, T)
    
    return M.*S    
end

function build_rhs_source(SD::NSD_2D,
                          QT::Exact,
                          q::Array,
                          mesh::St_mesh,
                          M::Matrix, #M is sparse for exact integration
                          T)

    S = user_source(q, mesh, T)
    
    return M*S   
end

