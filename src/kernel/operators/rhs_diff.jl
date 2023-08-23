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
# CompEuler
#
function build_rhs_diff_work_array()
    
    #qnel = zeros(mesh.ngl, mesh.nelem, neqs)
    ρel = zeros(mesh.ngl, mesh.nelem)
    uel = zeros(mesh.ngl, mesh.nelem)
    Tel = zeros(mesh.ngl, mesh.nelem)
    Eel = zeros(mesh.ngl, mesh.nelem)

    
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
        
        for j=1:mesh.ngl, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]
            
            ρel[i,j,iel] = qq[m,1]          
            uel[i,j,iel] = qq[m,2]/ρel[i,j,iel]
            vel[i,j,iel] = qq[m,3]/ρel[i,j,iel]
            
            Tel[i,j,iel] = qq[m,4]/ρel[i,j,iel] - δenergy*0.5*(uel[i,j,iel]^2 + vel[i,j,iel]^2)
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
    
end
