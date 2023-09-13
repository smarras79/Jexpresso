#---------------------------------------------------------------------------
# Optimized (more coud possibly be done)
#---------------------------------------------------------------------------
function build_rhs!(RHS, u, params, time)
    #
    # build_rhs()! is called by TimeIntegrators.jl -> time_loop!() via ODEProblem(rhs!, u, tspan, params)
    #
    _build_rhs!(RHS, u, params, time)
    
end

function RHStoDU!(du, RHS, neqs, npoin)
    for i=1:neqs
        idx = (i-1)*npoin
        du[idx+1:i*npoin] = @view RHS[:,i]
    end  
end

function u2uaux!(uaux, u, neqs, npoin)

    for i=1:neqs
        idx = (i-1)*npoin
        uaux[:,i] = view(u, idx+1:i*npoin)
    end
    
end


function uaux2u!(u, uaux, neqs, npoin)

    for i=1:neqs
        idx = (i-1)*npoin
        for j=1:npoin
            u[idx+j] = uaux[j,i]
        end
    end
    
end

function resetRHSToZero_inviscid!(params)
    fill!(params.rhs_el, zero(params.T))   
    fill!(params.RHS,    zero(params.T))
end

function resetRHSToZero_viscous!(params)
    fill!(params.rhs_diff_el,  zero(params.T))
    fill!(params.rhs_diffξ_el, zero(params.T))
    fill!(params.rhs_diffη_el, zero(params.T))
    fill!(params.RHS_visc,     zero(params.T))
end


function uToPrimitives!(uprimitive, u, mesh, δtotal_energy, iel)

    PhysConst = PhysicalConst{Float64}()

    for j=1:mesh.ngl, i=1:mesh.ngl
        
        m1 = mesh.connijk[iel,i,j]
        m2 = mesh.npoin + m1
        m3 = 2*mesh.npoin + m1
        m4 = 3*mesh.npoin + m1
        
        uprimitive[i,j,1] = u[m1]
        uprimitive[i,j,2] = u[m2]/u[m1]
        uprimitive[i,j,3] = u[m3]/u[m1]
        uprimitive[i,j,4] = u[m4]/u[m1] - δtotal_energy*0.5*(uprimitive[i,j,2]^2 + uprimitive[i,j,3]^2)

        #Pressure:
        uprimitive[i,j,end] = perfectGasLaw_ρθtoP(PhysConst, ρ=uprimitive[i,j,1], θ=uprimitive[i,j,4])
        
    end
    
end


function rhs!(du, u, params, time)
    
    build_rhs!(@view(params.RHS[:,:]), u, params, time)
    
    RHStoDU!(du, @view(params.RHS[:,:]), params.neqs, params.mesh.npoin)
end


function _build_rhs!(RHS, u, params, time)

    T       = Float64
    SD      = params.SD
    QT      = params.QT
    CL      = params.CL
    neqs    = params.neqs
    ngl     = params.mesh.ngl
    nelem   = params.mesh.nelem
    npoin   = params.mesh.npoin
    
    #-----------------------------------------------------------------------------------
    # Inviscid rhs:
    #-----------------------------------------------------------------------------------    
    resetRHSToZero_inviscid!(params) 
    
    inviscid_rhs_el!(u, params, true, SD)
    
    apply_boundary_conditions!(u, params.uaux, time,
                               params.mesh, params.metrics, params.basis,
                               params.rhs_el, params.ubdy,
                               params.ω, SD, neqs, params.inputs)
    
    DSS_rhs!(@view(params.RHS[:,:]), @view(params.rhs_el[:,:,:,:]), params.mesh, nelem, ngl, neqs, SD)
    
    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)

        resetRHSToZero_viscous!(params)
        
        viscous_rhs_el!(u, params, SD)
        
        DSS_rhs!(@view(params.RHS_visc[:,:]), @view(params.rhs_diff_el[:,:,:,:]), params.mesh, nelem, ngl, neqs, SD)
        
        params.RHS[:,:] .= @view(params.RHS[:,:]) .+ @view(params.RHS_visc[:,:])
    end
    
    divive_by_mass_matrix!(@view(params.RHS[:,:]), @view(params.M[:]), QT, neqs)
    
end


function inviscid_rhs_el!(u, params, lsource, SD::NSD_2D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    for iel=1:params.mesh.nelem
        
        uToPrimitives!(params.uprimitive, u, params.mesh, params.inputs[:δtotal_energy], iel)
        
        for j=1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,j]
            
            user_flux!(@view(params.F[i,j,:]), @view(params.G[i,j,:]), SD,
                       @view(params.uaux[ip,:]), 
                       params.qe[ip,end],         #pref
                       params.mesh; neqs=params.neqs)
            
            if lsource
                user_source!(@view(params.S[i,j,:]),
                             @view(params.uaux[ip,:]),
                             params.qe[ip,1],          #ρref 
                             params.mesh.npoin; neqs=params.neqs)
            end
        end

        #= WIP NONCON.
        _expansion_inviscid!(@view(params.rhs_el[iel,:,:,:]), params.uprimitive,
                             params.mesh, params.metrics, params.basis,
                             @view(params.F[:,:,:]), @view(params.G[:,:,:]), @view(params.S[:,:,:]),
                             params.ω, params.mesh.ngl, params.mesh.npoin, params.neqs, 0, iel,
                             params.CL, params.QT, SD)=#
        
        _expansion_inviscid!(params, iel, params.CL, params.QT, SD)
        
    end
end

function viscous_rhs_el!(u, params, SD::NSD_2D)
    
    for iel=1:params.mesh.nelem
        
        uToPrimitives!(params.uprimitive, u, params.mesh, params.inputs[:δtotal_energy], iel)

        for ieq=2:params.neqs
            _expansion_visc!(@view(params.rhs_diffξ_el[iel,:,:,ieq]), @view(params.rhs_diffη_el[iel,:,:,ieq]), @view(params.uprimitive[:,:,ieq]), params.visc_coeff[ieq], params.ω, params.mesh, params.basis, params.metrics, params.inputs, iel, ieq, params.QT, SD)
        end
        
    end
    params.rhs_diff_el .= @views (params.rhs_diffξ_el .+ params.rhs_diffη_el)
end


function _expansion_inviscid!(rhs_el, uprimitive, mesh, metrics, basis, F, G, S, ω, ngl, npoin, neqs, ieq, iel, ::NCL, QT::Inexact, SD::NSD_2D)
    
    for ieq=1:neqs
        for j=1:ngl
            for i=1:ngl
                ωJac = ω[i]*ω[j]*metrics.Je[iel,i,j]
                
                dFdξ = 0.0; dFdη = 0.0
                dGdξ = 0.0; dGdη = 0.0
                dpdξ = 0.0; dpdη = 0.0               
                for k = 1:ngl
                    dFdξ += basis.dψ[k,i]*F[k,j,ieq]
                    dFdη += basis.dψ[k,j]*F[i,k,ieq]
                    
                    dGdξ += basis.dψ[k,i]*G[k,j,ieq]
                    dGdη += basis.dψ[k,j]*G[i,k,ieq]
                                        
                    dpdξ += basis.dψ[k,i]*uprimitive[k,j,neqs+1]
                    dpdη += basis.dψ[k,j]*uprimitive[i,k,neqs+1]
                end
                dξdx_ij = metrics.dξdx[iel,i,j]
                dξdy_ij = metrics.dξdy[iel,i,j]
                dηdx_ij = metrics.dηdx[iel,i,j]
                dηdy_ij = metrics.dηdy[iel,i,j]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij            
                dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij

                dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij            
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
                dpdx = dpdξ*dξdx_ij + dpdη*dηdx_ij            
                dpdy = dpdξ*dξdy_ij + dpdη*dηdy_ij

                ρij = uprimitive[i,j,1]
                uij = uprimitive[i,j,2]
                vij = uprimitive[i,j,3]
                
                if (ieq == 1)
                    auxi = ωJac*(dFdx + dGdy)
                elseif(ieq == 2)
                    auxi = ωJac*(uij*dFdx + vij*dGdy + dpdx/ρij)
                elseif(ieq == 3)
                    auxi = ωJac*(uij*dFdx + vij*dGdy + dpdy/ρij )#- S[i,j,ieq])
                elseif(ieq == 4)
                    auxi = ωJac*(uij*dFdx + vij*dGdy)
                end
                
                rhs_el[i,j,ieq] -= auxi
            end
        end
    end        
end

function _expansion_inviscid!(params, iel, ::CL, QT::Inexact, SD::NSD_2D)

    for ieq=1:params.neqs
        for j=1:params.mesh.ngl
            for i=1:params.mesh.ngl
                ωJac = params.ω[i]*params.ω[j]*params.metrics.Je[iel,i,j]
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                @turbo for k = 1:params.mesh.ngl
                    dFdξ += params.basis.dψ[k,i]*params.F[k,j,ieq]
                    dFdη += params.basis.dψ[k,j]*params.F[i,k,ieq]
                    
                    dGdξ += params.basis.dψ[k,i]*params.G[k,j,ieq]
                    dGdη += params.basis.dψ[k,j]*params.G[i,k,ieq]
                end
                dξdx_ij = params.metrics.dξdx[iel,i,j]
                dξdy_ij = params.metrics.dξdy[iel,i,j]
                dηdx_ij = params.metrics.dηdx[iel,i,j]
                dηdy_ij = params.metrics.dηdy[iel,i,j]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
                dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij

                dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
                auxi = ωJac*((dFdx + dGdy) - params.S[i,j,ieq])
                params.rhs_el[iel,i,j,ieq] -= auxi
            end
        end
    end
end


function _expansion_inviscid!(rhs_el, uprimitiveieq, mesh, metrics, basis, F, G, S, ω, ngl, npoin, neqs, ieq, iel, ::CL, QT::Exact, SD::NSD_2D)
    
    N = ngl
    Q = N + 1
    for l=1:Q
        for k=1:Q
            ωJac = ω[k]*ω[l]*metrics.Je[iel,k,l]
            
            dFdξ = 0.0
            dFdη = 0.0
            dGdξ = 0.0
            dGdη = 0.0
            for n = 1:N
                for m = 1:N
                    dFdξ += basis.dψ[m,k]* basis.ψ[n,l]*F[m,n]
                    dFdη +=  basis.ψ[m,k]*basis.dψ[n,l]*F[m,n]
                    
                    dGdξ += basis.dψ[m,k]* basis.ψ[n,l]*G[m,n]
                    dGdη +=  basis.ψ[m,k]*basis.dψ[n,l]*G[m,n]
                end
            end
            
            dξdx_kl = metrics.dξdx[iel,k,l]
            dξdy_kl = metrics.dξdy[iel,k,l]
            dηdx_kl = metrics.dηdx[iel,k,l]
            dηdy_kl = metrics.dηdy[iel,k,l]
            for j = 1:N
                for i = 1:N
                    #I = mesh.connijk[iel,i,j]
                    
                    dFdx = dFdξ*dξdx_kl + dFdη*dηdx_kl
                    dGdx = dGdξ*dξdx_kl + dGdη*dηdx_kl

                    dFdy = dFdξ*dξdy_kl + dFdη*dηdy_kl
                    dGdy = dGdξ*dξdy_kl + dGdη*dηdy_kl
                    
                    #auxi = ωJac*((dFdx + dGdy) - S[i,j])
                    auxi = ωJac*basis.ψ[i,k]*basis.ψ[j,l]*((dFdx + dGdy) - S[i,j])
                    rhs_el[i,j] -= auxi
                    #RHS[I] -= auxi
                end
            end
        end
    end    
end


function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeffieq, ω, mesh, basis, metrics, inputs, iel, ieq, QT::Inexact, SD::NSD_2D)
  
    for l = 1:mesh.ngl
        for k = 1:mesh.ngl
            ωJac = ω[k]*ω[l]*metrics.Je[iel,k,l]
            
            dudξ = 0.0
            dudη = 0.0
            @turbo for ii = 1:mesh.ngl
                dudξ += basis.dψ[ii,k]*uprimitiveieq[ii,l]
                dudη += basis.dψ[ii,l]*uprimitiveieq[k,ii]
            end
            dξdx_kl = metrics.dξdx[iel,k,l]
            dξdy_kl = metrics.dξdy[iel,k,l]
            dηdx_kl = metrics.dηdx[iel,k,l]
            dηdy_kl = metrics.dηdy[iel,k,l]
            
            auxi = dudξ*dξdx_kl + dudη*dηdx_kl
            dudx = visc_coeffieq*auxi
            
            auxi = dudξ*dξdy_kl + dudη*dηdy_kl
            dudy = visc_coeffieq*auxi
            
            ∇ξ∇u_kl = (dξdx_kl*dudx + dξdy_kl*dudy)*ωJac
            ∇η∇u_kl = (dηdx_kl*dudx + dηdy_kl*dudy)*ωJac     
            
            @turbo for i = 1:mesh.ngl
                dhdξ_ik = basis.dψ[i,k]
                dhdη_il = basis.dψ[i,l]
                
                rhs_diffξ_el[i,l] -= dhdξ_ik * ∇ξ∇u_kl
                rhs_diffη_el[k,i] -= dhdη_il * ∇η∇u_kl
            end
        end
    end  
end
