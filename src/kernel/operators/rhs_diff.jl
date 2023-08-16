#---------------------------------------------------------------------------
# Fetch equations name to access the user_rhs functions
#---------------------------------------------------------------------------
if (length(ARGS) === 1) #equations
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "/../../equations/", ARGS[1], "/user_source.jl")
    else
        user_source_dir = "../../fallbacks/source.jl"
    end
elseif (length(ARGS) === 2)  #equations/equations_case_name
    user_flux_dir   = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_flux.jl")
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl")
    else
        @info " user_source.jl not defined. The fallback ../../fallbacks/source.jl will be used."
        user_source_dir =  "../../fallbacks/source.jl"
    end
end
include(user_flux_dir)
include(user_source_dir)
include("../ArtificialViscosity/DynSGS.jl")
#---------------------------------------------------------------------------
#
# AdvDiff
#
function build_rhs_diff(SD::NSD_1D, QT::Inexact, PT::AdvDiff, qp::Array, nvars, basis, ω, inputs,  mesh::St_mesh, metrics::St_metrics, μ, T;)

    N           = mesh.ngl - 1
    qnel        = zeros(mesh.ngl, mesh.nelem)
    rhsdiffξ_el = zeros(mesh.ngl, mesh.nelem)
    
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2.0
        dξdx = 2.0/mesh.Δx[iel]
        
        for i=1:mesh.ngl
            qnel[i,iel,1] = qp[mesh.connijk[i,iel], 1]
        end
        
        for k = 1:mesh.ngl
            ωJk = ω[k]*Jac
            
            dqdξ = 0.0
            for i = 1:mesh.ngl
                dqdξ = dqdξ + basis.dψ[i,k]*qnel[i,iel]
            end
            dqdx = dqdξ*dξdx            
            ∇ξ∇q = dξdx*dqdx
            
            for i = 1:mesh.ngl
                hll     = basis.ψ[k,k]
                dhdξ_ik = basis.dψ[i,k]
                
                rhsdiffξ_el[i, iel] -= ωJk * basis.dψ[i,k] * basis.ψ[k,k]*∇ξ∇q
            end
        end
    end
      
    return rhsdiffξ_el*inputs[:νx]
end

function build_rhs_diff(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, nvars, basis, ω, inputs,  mesh::St_mesh, metrics::St_metrics, μ, T;)

    N = mesh.ngl - 1
    
    qnel = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    
    rhsdiffξ_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    rhsdiffη_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]            
            qnel[i,j,iel] = qp[m,1]
        end
        
        for k = 1:mesh.ngl, l = 1:mesh.ngl
            ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]
            
            dqdξ = 0.0
            dqdη = 0.0
            for i = 1:mesh.ngl
                dqdξ = dqdξ + basis.dψ[i,k]*qnel[i,l,iel]
                dqdη = dqdη + basis.dψ[i,l]*qnel[k,i,iel]
            end
            
            dqdx = dqdξ*metrics.dξdx[k,l,iel] + dqdη*metrics.dηdx[k,l,iel]
            dqdy = dqdξ*metrics.dξdy[k,l,iel] + dqdη*metrics.dηdy[k,l,iel]
            
            ∇ξ∇q_kl = metrics.dξdx[k,l,iel]*dqdx + metrics.dξdy[k,l,iel]*dqdy
            ∇η∇q_kl = metrics.dηdx[k,l,iel]*dqdx + metrics.dηdy[k,l,iel]*dqdy
            
            for i = 1:mesh.ngl
                hll,     hkk     =  basis.ψ[l,l],  basis.ψ[k,k]
                dhdξ_ik, dhdη_il = basis.dψ[i,k], basis.dψ[i,l]
                
                rhsdiffξ_el[i,l,iel] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
                rhsdiffη_el[k,i,iel] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
            end
        end
    end

    return (rhsdiffξ_el*inputs[:νx] + rhsdiffη_el*inputs[:νy])
    
end

#
# LinearCLaw
#
function build_rhs_diff(SD::NSD_2D, QT, PT::LinearCLaw, qp, neqs, basis, ω, inputs,  mesh::St_mesh, metrics::St_metrics, μ, T;)
    
    N = mesh.ngl - 1

    qnel = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)

    rhsdiffξ_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqs)
    rhsdiffη_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqs)
    qq = zeros(mesh.npoin,neqs)

    #
    # qp[1:npoin]         <-- qq[1:npoin, "p"]
    # qp[npoin+1:2npoin]  <-- qq[1:npoin, "u"]
    # qp[2npoin+1:3npoin] <-- qq[1:npoin, "v"]
    #
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] = qp[idx+1:i*mesh.npoin]
    end
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]
            qnel[i,j,iel,1:neqs] = qq[m,1:neqs]
        end

        for k = 1:mesh.ngl, l = 1:mesh.ngl
            ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]

            for ieq = 1:neqs
                dqdξ = 0.0
                dqdη = 0.0
                for i = 1:mesh.ngl
                    dqdξ = dqdξ + basis.dψ[i,k]*qnel[i,l,iel,ieq]
                    dqdη = dqdη + basis.dψ[i,l]*qnel[k,i,iel,ieq]
                end
                dqdx = dqdξ*metrics.dξdx[k,l,iel] + dqdη*metrics.dηdx[k,l,iel]
                dqdy = dqdξ*metrics.dξdy[k,l,iel] + dqdη*metrics.dηdy[k,l,iel]

                ∇ξ∇q_kl = metrics.dξdx[k,l,iel]*dqdx + metrics.dξdy[k,l,iel]*dqdy
                ∇η∇q_kl = metrics.dηdx[k,l,iel]*dqdx + metrics.dηdy[k,l,iel]*dqdy

                for i = 1:mesh.ngl

                    hll,     hkk     = basis.ψ[l,l],  basis.ψ[k,k]
                    dhdξ_ik, dhdη_il = basis.dψ[i,k], basis.dψ[i,l]

                    rhsdiffξ_el[i,l,iel, ieq] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
                    rhsdiffη_el[k,i,iel, ieq] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
                end
            end
        end
     end

    return (rhsdiffξ_el*inputs[:νx] + rhsdiffη_el*inputs[:νy])

end


#
# ShallowWater
#
function build_rhs_diff(SD::NSD_1D, QT, PT::ShallowWater, qp, neqs, basis, ω, inputs,  mesh::St_mesh, metrics::St_metrics, μ, T;)

    N = mesh.ngl - 1

    qnel = zeros(mesh.ngl, mesh.nelem, neqs)

    rhsdiffξ_el = zeros(mesh.ngl, mesh.nelem, neqs)
    qq = zeros(mesh.npoin,neqs)

    #
    # qp[1:npoin]         <-- qq[1:npoin, "p"]
    # qp[npoin+1:2npoin]  <-- qq[1:npoin, "u"]
    # qp[2npoin+1:3npoin] <-- qq[1:npoin, "v"]
    #
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] = qp[idx+1:i*mesh.npoin]
    end
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2.0
        for i=1:mesh.ngl
            m = mesh.conn[i,iel]
            qnel[i,iel,1] = qq[m,1]
            qnel[i,iel,2] = qq[m,2]/qq[m,1]
        end
        dξdx = 2.0/mesh.Δx[iel]
        for k = 1:mesh.ngl
            ωJkl = ω[k]*Jac

            for ieq = 1:neqs
                dqdξ = 0.0
                for i = 1:mesh.ngl
                    dqdξ = dqdξ + basis.dψ[i,k]*qnel[i,iel,ieq]
                    #@info "contribution", basis.dψ[i,k]*qnel[k,iel,ieq]
                end
                #@info "dqdxi", dqdξ
                #if (ieq > 1)
                dqdx = μ[iel] * (dqdξ) * dξdx
                #else
                #    dqdx = 0.0
                #end 
                #@info "dqdx", dqdx, "vx", νx
                
                if (ieq > 1)
                    ip = mesh.conn[k,iel]
                    x = mesh.x[ip]
                    Hb = bathymetry(x)
                    Hs = max(qq[ip,1] - Hb,0.001)
                    dqdx = dqdx * qq[ip,1]#* Hs
                end

                ∇ξ∇q_kl =  dqdx*dξdx 
                for i = 1:mesh.ngl

                    hkk     = basis.ψ[k,k]
                    dhdξ_ik = basis.dψ[i,k]

                    rhsdiffξ_el[i,iel,ieq] -= ωJkl*dhdξ_ik*hkk*∇ξ∇q_kl
                end
            end
        end
     end
    
    return (rhsdiffξ_el)

end

function build_rhs_diff(SD::NSD_2D, QT, PT::ShallowWater, qp, neqs, basis, ω, inputs,  mesh::St_mesh, metrics::St_metrics, T;)
    
    N = mesh.ngl - 1

    qnel = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)

    rhsdiffξ_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqs)
    rhsdiffη_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqs)
    qq = zeros(mesh.npoin,neqs)

    #
    # qp[1:npoin]         <-- qq[1:npoin, "p"]
    # qp[npoin+1:2npoin]  <-- qq[1:npoin, "u"]
    # qp[2npoin+1:3npoin] <-- qq[1:npoin, "v"]
    #
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] = qp[idx+1:i*mesh.npoin]
    end
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]
            qnel[i,j,iel,1] = qq[m,1]
            qnel[i,j,iel,2] = qq[m,2]/qq[m,1]
            qnel[i,j,iel,3] = qq[m,3]/qq[m,1]
        end

        for k = 1:mesh.ngl, l = 1:mesh.ngl
            ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]

            for ieq = 1:neqs
                dqdξ = 0.0
                dqdη = 0.0
                for i = 1:mesh.ngl
                    dqdξ = dqdξ + basis.dψ[i,k]*qnel[i,l,iel,ieq]
                    dqdη = dqdη + basis.dψ[i,l]*qnel[k,i,iel,ieq]
                end
                dqdx = inputs[:νx] * (dqdξ*metrics.dξdx[k,l,iel] + dqdη*metrics.dηdx[k,l,iel])
                dqdy = inputs[:νy] * (dqdξ*metrics.dξdy[k,l,iel] + dqdη*metrics.dηdy[k,l,iel])
                if (ieq > 1)
                    ip = mesh.connijk[k,l,iel]
                    x = mesh.x[ip]
                    y = mesh.y[ip]
                    Hb = bathymetry(x,y)
                    Hs = qq[ip,1] - Hb
                    dqdx = dqdx * Hs
                    dqdy = dqdy * Hs
                end                

                ∇ξ∇q_kl =  (metrics.dξdx[k,l,iel]*dqdx + metrics.dξdy[k,l,iel]*dqdy)
                ∇η∇q_kl =  (metrics.dηdx[k,l,iel]*dqdx + metrics.dηdy[k,l,iel]*dqdy)

                for i = 1:mesh.ngl

                    hll,     hkk     = basis.ψ[l,l],  basis.ψ[k,k]
                    dhdξ_ik, dhdη_il = basis.dψ[i,k], basis.dψ[i,l]
     
                    rhsdiffξ_el[i,l,iel, ieq] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
                    rhsdiffη_el[k,i,iel, ieq] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
                end
            end
        end
     end

    return (rhsdiffξ_el + rhsdiffη_el)

end

#
# CompEuler
#
function build_rhs_diff_work_array()
    
    #qnel = zeros(mesh.ngl, mesh.nelem, neqs)
    ρel = zeros(mesh.ngl, mesh.nelem)
    uel = zeros(mesh.ngl, mesh.nelem)
    Tel = zeros(mesh.ngl, mesh.nelem)
    Eel = zeros(mesh.ngl, mesh.nelem)

    
end

function build_rhs_diff(SD::NSD_1D, QT, PT::CompEuler, qp, neqs, basis, ω, inputs,  mesh::St_mesh, metrics::St_metrics, μ, T;)
    
    ρel = zeros(T, mesh.ngl, mesh.nelem)
    uel = zeros(T, mesh.ngl, mesh.nelem)
        
    rhsdiffξ_el = zeros(T, mesh.ngl, mesh.nelem, neqs)
   
    PhysConst = PhysicalConst{Float64}()
    γ  = PhysConst.γ
    Pr = PhysConst.Prnum
    
    qq = zeros(T, mesh.npoin,neqs)
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= 0.0 .+ view(qp, idx+1:i*mesh.npoin)
    end
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2.0
        for i=1:mesh.ngl
            m = mesh.conn[i,iel]
             
            ρel[i,iel] = qq[m,1]
        end
        
        ν = Pr*μ[iel]/maximum(ρel[:,iel])
        κ = Pr*μ[iel]/(γ - 1.0)
        
        dξdx = 2.0/mesh.Δx[iel]
        for k = 1:mesh.ngl
            ωJkl = ω[k]*Jac

            dρdξ = 0.0
            dudξ = 0.0
            dTdξ = 0.0
            for i = 1:mesh.ngl

                m = mesh.conn[i,iel]
                uel[i,iel] = qq[m,2]/qq[m,1]
                Tel        = qq[m,3]/qq[m,1] - 0.5*uel[i,iel]^2
                
                dρdξ = dρdξ + basis.dψ[i,k]*ρel[i,iel]
                dudξ = dudξ + basis.dψ[i,k]*uel[i,iel]
                dTdξ = dTdξ + basis.dψ[i,k]*Tel
            end
             
            dρdx =  ν * dρdξ*dξdx
            dudx =  μ[iel] * dudξ*dξdx
            dTdx =  μ[iel] * dudξ*dξdx * uel[k,iel] + κ * dTdξ*dξdx
            
            ∇ξ∇ρ_kl = dρdx*dξdx
            ∇ξ∇u_kl = dudx*dξdx
            ∇ξ∇T_kl = dTdx*dξdx
            for i = 1:mesh.ngl
                
                rhsdiffξ_el[i,iel,1] -= ωJkl*basis.dψ[i,k]*∇ξ∇ρ_kl
                rhsdiffξ_el[i,iel,2] -= ωJkl*basis.dψ[i,k]*∇ξ∇u_kl
                rhsdiffξ_el[i,iel,3] -= ωJkl*basis.dψ[i,k]*∇ξ∇T_kl
            end
        end
    end
    
    return rhsdiffξ_el

end

function flux2primitives(q)
    
end

function build_rhs_diff!(rhs_diff::SubArray{Float64}, SD::NSD_2D, QT, PT::CompEuler, qp, neqs, basis, ω, inputs, mesh::St_mesh, metrics::St_metrics, μ, T; qoutauxi=zeros(1,1))
    
    ρel = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    uel = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    vel = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    Tel = zeros(mesh.ngl, mesh.ngl, mesh.nelem)

    rhsdiffξ_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqs)
    rhsdiffη_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqs)
    
    qq = zeros(mesh.npoin, neqs)
    
    #
    # qp[1:npoin]         <-- qq[1:npoin, "ρ"]
    # qp[npoin+1:2npoin]  <-- qq[1:npoin, "ρu"]
    # qp[2npoin+1:3npoin] <-- qq[1:npoin, "ρE"]
    #
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] = qp[idx+1:i*mesh.npoin]
    end
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    if (inputs[:case] === "rtb")
        δenergy = 0.0
    else
        δenergy = 1.0
    end

    
    γ = 1.4
    Pr = 0.1
    for iel=1:mesh.nelem
        
        @inbounds begin
            for j=1:mesh.ngl, i=1:mesh.ngl
                m = mesh.connijk[i,j,iel]
                
                ρel[i,j,iel] = qq[m,1]          
                uel[i,j,iel] = qq[m,2]/ρel[i,j,iel]
                vel[i,j,iel] = qq[m,3]/ρel[i,j,iel]
                
                Tel[i,j,iel] = qq[m,4]/ρel[i,j,iel] - δenergy*0.5*(uel[i,j,iel]^2 + vel[i,j,iel]^2)
            end
        end
        #ν = Pr*μ[iel]/maximum(ρel[:,:,iel])
        #κ = Pr*μ[iel]/(γ - 1.0)
        #ν = μ[iel]#10.0
        #κ = μ[iel]#10.0
        
        ν = 0.0
        μ[iel] = inputs[:νx]
        κ = μ[iel]
        for l = 1:mesh.ngl
            for k = 1:mesh.ngl
                ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]
                
                dρdξ = 0.0
                dudξ = 0.0
                dvdξ = 0.0
                dTdξ = 0.0

                dρdη = 0.0
                dudη = 0.0
                dvdη = 0.0
                dTdη = 0.0
                for i = 1:mesh.ngl
                    dρdξ += basis.dψ[i,k]*ρel[i,l,iel]
                    dudξ += basis.dψ[i,k]*uel[i,l,iel]
                    dvdξ += basis.dψ[i,k]*vel[i,l,iel]
                    dTdξ += basis.dψ[i,k]*Tel[i,l,iel]

                    dρdη += basis.dψ[i,l]*ρel[k,i,iel]
                    dudη += basis.dψ[i,l]*uel[k,i,iel]
                    dvdη += basis.dψ[i,l]*vel[k,i,iel]
                    dTdη += basis.dψ[i,l]*Tel[k,i,iel]
                end
                                
                #
                dξdx_kl = metrics.dξdx[k,l,iel]
                dξdy_kl = metrics.dξdy[k,l,iel]
                dηdx_kl = metrics.dηdx[k,l,iel]
                dηdy_kl = metrics.dηdy[k,l,iel]
                
                #
                dρdx =       ν*(dρdξ*dξdx_kl + dρdη*dηdx_kl)
                dudx =  μ[iel]*(dudξ*dξdx_kl + dudη*dηdx_kl)
                dvdx =  μ[iel]*(dvdξ*dξdx_kl + dvdη*dηdx_kl)
                dTdx =       κ*(dTdξ*dξdx_kl + dTdη*dηdx_kl) #+μ∇u⋅u
                
                dρdy =       ν*(dρdξ*dξdy_kl + dρdη*dηdy_kl)
                dudy =  μ[iel]*(dudξ*dξdy_kl + dudη*dηdy_kl)
                dvdy =  μ[iel]*(dvdξ*dξdy_kl + dvdη*dηdy_kl)
                dTdy =       κ*(dTdξ*dξdy_kl + dTdη*dηdy_kl) #+μ∇u⋅u
                
                ∇ξ∇ρ_kl = dξdx_kl*dρdx + dξdy_kl*dρdy
                ∇η∇ρ_kl = dηdx_kl*dρdx + dηdy_kl*dρdy
                
                ∇ξ∇u_kl = dξdx_kl*dudx + dξdy_kl*dudy
                ∇η∇u_kl = dηdx_kl*dudx + dηdy_kl*dudy            
                ∇ξ∇v_kl = dξdx_kl*dvdx + dξdy_kl*dvdy
                ∇η∇v_kl = dηdx_kl*dvdx + dηdy_kl*dvdy

                ∇ξ∇T_kl = dξdx_kl*dTdx + dξdy_kl*dTdy
                ∇η∇T_kl = dηdx_kl*dTdx + dηdy_kl*dTdy
                
                for i = 1:mesh.ngl
                    
                    dhdξ_ik, dhdη_il = basis.dψ[i,k], basis.dψ[i,l]
                    
                    rhsdiffξ_el[i,l,iel,1] -= ωJkl*dhdξ_ik*∇ξ∇ρ_kl
                    rhsdiffη_el[k,i,iel,1] -= ωJkl*dhdη_il*∇η∇ρ_kl
                    
                    rhsdiffξ_el[i,l,iel,2] -= ωJkl*dhdξ_ik*∇ξ∇u_kl
                    rhsdiffη_el[k,i,iel,2] -= ωJkl*dhdη_il*∇η∇u_kl
                    
                    rhsdiffξ_el[i,l,iel,3] -= ωJkl*dhdξ_ik*∇ξ∇v_kl
                    rhsdiffη_el[k,i,iel,3] -= ωJkl*dhdη_il*∇η∇v_kl
                    
                    rhsdiffξ_el[i,l,iel,4] -= ωJkl*dhdξ_ik*∇ξ∇T_kl
                    rhsdiffη_el[k,i,iel,4] -= ωJkl*dhdη_il*∇η∇T_kl
                    
                end
            end
        end
    end

    rhs_diff .= @views (rhsdiffξ_el[:,:,:,:] + rhsdiffη_el[:,:,:,:])
#return (rhsdiffξ_el + rhsdiffη_el)

end
