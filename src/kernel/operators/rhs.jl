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
    fill!(params.b,      zero(params.T))
    fill!(params.B,      zero(params.T))
end

function reset_laguerre_filters!(params)
    fill!(params.b_lag,      zero(params.T))
    fill!(params.B_lag,      zero(params.T))
end

function resetRHSToZero_viscous!(params)
    fill!(params.rhs_diff_el,  zero(params.T))
    fill!(params.rhs_diffξ_el, zero(params.T))
    fill!(params.rhs_diffη_el, zero(params.T))
    fill!(params.RHS_visc,     zero(params.T))
end

function uToPrimitives!(neqs, uprimitive, u, uauxe, mesh, δtotal_energy, iel, PT, ::CL, ::AbstractPert, SD::NSD_1D)
    
    γ = 1.4
    γm1 = 0.4
    for i=1:mesh.ngl
        
        m1 = mesh.connijk[iel,i,1]
        m2 = m1 + mesh.npoin
        m3 = m2 + mesh.npoin
        
        x  = mesh.x[m1]           
        A = 1.0 + 2.2*(x - 1.5)^2
        
        uprimitive[i,1,1] = u[m1]/A #ρ
        uprimitive[i,1,2] = u[m2]/u[m1] #u
        uprimitive[i,1,3] = γm1*(u[m3]/u[m1] - 0.5*γ*uprimitive[i,j,1]*uprimitive[i,j,1]) #T

        #Pressure:
        uprimitive[i,1,end] = uprimitive[i,1,1]*uprimitive[i,1,3] #ρ*T
           
    end
end

function uToPrimitives!(neqs, uprimitive, u, uauxe, mesh, δtotal_energy, iel, PT, ::CL, ::TOTAL, SD::NSD_2D)

    if typeof(PT) == CompEuler
        
        PhysConst = PhysicalConst{Float64}()
        
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            m1 = mesh.connijk[iel,i,j]
            m2 = m1 + mesh.npoin
            m3 = m2 + mesh.npoin
            m4 = m3 + mesh.npoin
            
            uprimitive[i,j,1] = u[m1]
            uprimitive[i,j,2] = u[m2]/u[m1]
            uprimitive[i,j,3] = u[m3]/u[m1]
            uprimitive[i,j,4] = u[m4]/u[m1] - δtotal_energy*0.5*(uprimitive[i,j,2]^2 + uprimitive[i,j,3]^2)
            #Tracers
            if(neqs > 4)
                mieq = m4
                for ieq = 5:neqs
                    mieq = mieq + mesh.npoin
                    uprimitive[i,j,ieq] = u[mieq]
                end
            end
            #Pressure:
            uprimitive[i,j,end] = perfectGasLaw_ρθtoP(PhysConst, ρ=uprimitive[i,j,1], θ=uprimitive[i,j,4])
            
        end
        
    elseif typeof(PT) == AdvDiff
        
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            mieq = mesh.connijk[iel,i,j]
            uprimitive[i,j,1] = u[mieq]
            
            for ieq = 2:neqs
                mieq = mieq + mesh.npoin
                uprimitive[i,j,ieq] = u[mieq]
            end
            
        end
    end
end

function uToPrimitives!(neqs, uprimitive, u, uauxe, mesh, δtotal_energy, iel, PT, ::CL, ::PERT, SD::NSD_2D)
    
    if typeof(PT) == CompEuler
        PhysConst = PhysicalConst{Float64}()
        
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            m1 = mesh.connijk[iel,i,j]
            m2 = m1 + mesh.npoin
            m3 = m2 + mesh.npoin
            m4 = m3 + mesh.npoin

            uprimitive[i,j,1] = u[m1] + uauxe[m1,1]
            uprimitive[i,j,2] = u[m2]/uprimitive[i,j,1]
            uprimitive[i,j,3] = u[m3]/uprimitive[i,j,1]
            uprimitive[i,j,4] = (u[m4] + uprimitive[i,j,1]*uauxe[m1,4])/uprimitive[i,j,1] # CHECK THIS FOR ENE- δtotal_energy*0.5*(uprimitive[i,j,2]^2 + uprimitive[i,j,3]^2)

            #Tracers
            if(neqs > 4)
                mieq = m4
                for ieq = 5:neqs
                    mieq = mieq + mesh.npoin
                    uprimitive[i,j,ieq] = u[mieq] + uauxe[mieq,1]
                end
            end
            
            #Pressure:
            uprimitive[i,j,end] = perfectGasLaw_ρθtoP(PhysConst, ρ=uprimitive[i,j,1], θ=uprimitive[i,j,4])
            
        end
        
    elseif typeof(PT) == AdvDiff
        
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            mieq = mesh.connijk[iel,i,j]
            uprimitive[i,j,1] = u[mieq]
            
            for ieq = 2:neqs
                mieq = mieq + mesh.npoin
                uprimitive[i,j,ieq] = u[mieq]
            end
            
        end
    end
    
end


function uToPrimitives!(neqs, uprimitive, u, uprimitivee, mesh, δtotal_energy, iel, PT, ::NCL, ::AbstractPert, SD::NSD_2D)

    if typeof(PT) == CompEuler
        PhysConst = PhysicalConst{Float64}()
        
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            m1 = mesh.connijk[iel,i,j]
            m2 = m1 + mesh.npoin
            m3 = m2 + mesh.npoin
            m4 = m3 + mesh.npoin

            uprimitive[i,j,1] = u[m1]
            uprimitive[i,j,2] = u[m2]
            uprimitive[i,j,3] = u[m3]
            uprimitive[i,j,4] = u[m4]

            #Tracers
            if(neqs > 4)
                mieq = 4
                for ieq = 5:neqs
                    mieq = mieq + mesh.npoin
                    uprimitive[i,j,ieq] = u[mieq]                
                end
            end
            
            #Pressure:
            uprimitive[i,j,end] = perfectGasLaw_ρθtoP(PhysConst, ρ=uprimitive[i,j,1], θ=uprimitive[i,j,4])
            
        end
        
    elseif typeof(PT) == AdvDiff
        
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            mieq = mesh.connijk[iel,i,j]
            uprimitive[i,j,1] = u[mieq]
            
            for ieq = 2:neqs
                mieq = mieq + mesh.npoin
                uprimitive[i,j,ieq] = u[mieq]
            end
            
        end
    end
    
end


function rhs!(du, u, params, time)
    
    build_rhs!(@view(params.RHS[:,:]), u, params, time)
    if (params.laguerre) 
        build_rhs_laguerre!(@view(params.RHS_lag[:,:]), u, params, time)
        params.RHS .= @views(params.RHS .+ params.RHS_lag)
    end
    RHStoDU!(du, @view(params.RHS[:,:]), params.neqs, params.mesh.npoin)
    
end

function _build_rhs!(RHS, u, params, time)

    T       = Float64
    SD      = params.SD
    QT      = params.QT
    CL      = params.CL
    AD      = params.AD
    neqs    = params.neqs
    ngl     = params.mesh.ngl
    nelem   = params.mesh.nelem
    npoin   = params.mesh.npoin
    lsource = params.inputs[:lsource]
        
    #-----------------------------------------------------------------------------------
    # Inviscid rhs:
    #-----------------------------------------------------------------------------------    
    resetRHSToZero_inviscid!(params) 
    if (params.laguerre && params.inputs[:lfilter])
         reset_laguerre_filters!(params)
    end
    if (params.inputs[:lfilter])
        filter!(u, params, time, SD,params.SOL_VARS_TYPE)
    end
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    apply_boundary_conditions!(u, params.uaux, time, params.qe,
                               params.mesh, params.metrics, params.basis,
                               params.RHS, params.rhs_el, params.ubdy,
                               params.ω, neqs, params.inputs, SD)
    
    inviscid_rhs_el!(u, params, lsource, SD)
    DSS_rhs!(@view(params.RHS[:,:]), @view(params.rhs_el[:,:,:,:]), params.mesh, nelem, ngl, neqs, SD, AD)
    
    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)
        
        resetRHSToZero_viscous!(params)
        
        viscous_rhs_el!(u, params, SD)
        
        DSS_rhs!(@view(params.RHS_visc[:,:]), @view(params.rhs_diff_el[:,:,:,:]), params.mesh, nelem, ngl, neqs, SD, AD)
        
        params.RHS[:,:] .= @view(params.RHS[:,:]) .+ @view(params.RHS_visc[:,:])
    end
    for ieq=1:neqs
        divide_by_mass_matrix!(@view(params.RHS[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
    end
    
end

function inviscid_rhs_el!(u, params, lsource, SD::NSD_1D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    xmax = params.xmax
    xmin = params.xmin
    ymax = params.ymax    
    for iel=1:params.mesh.nelem

        uToPrimitives!(params.neqs, params.uprimitive, u, params.qe, params.mesh, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)
        
        for i=1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,1]
            
            user_flux!(@view(params.F[i,1,:]), @view(params.G[i,1,:]), SD,
                       @view(params.uaux[ip,:]),
                       @view(params.qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs)
            
            if lsource
                user_source!(@view(params.S[i,1,:]),
                             @view(params.uaux[ip,:]),
                             @view(params.qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE; neqs=params.neqs, x=params.mesh.x[ip],y=params.mesh.y[ip],xmax=xmax,xmin=xmin,ymax=ymax)
            end
        end
        
        _expansion_inviscid!(params, iel, params.CL, params.QT, SD, params.AD)
        
    end
end

function inviscid_rhs_el!(u, params, lsource, SD::NSD_2D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    xmax = params.xmax
    xmin = params.xmin
    ymax = params.ymax    
    for iel = 1:params.mesh.nelem

        uToPrimitives!(params.neqs, params.uprimitive, u, params.qe, params.mesh, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)

        for j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,j]
            
            user_flux!(@view(params.F[i,j,:]), @view(params.G[i,j,:]), SD,
                       @view(params.uaux[ip,:]),
                       @view(params.qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs)
            
            if lsource
                user_source!(@view(params.S[i,j,:]),
                             @view(params.uaux[ip,:]),
                             @view(params.qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE; neqs=params.neqs, x=params.mesh.x[ip],y=params.mesh.y[ip],xmax=xmax,xmin=xmin,ymax=ymax)
            end
        end

        _expansion_inviscid!(params, iel, params.CL, params.QT, SD, params.AD)
        
    end
end

function viscous_rhs_el!(u, params, SD::NSD_2D)
    
    for iel=1:params.mesh.nelem
        
        uToPrimitives!(params.neqs, params.uprimitive, u, params.qe, params.mesh, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)

        for ieq in params.ivisc_equations
            _expansion_visc!(@view(params.rhs_diffξ_el[iel,:,:,ieq]), @view(params.rhs_diffη_el[iel,:,:,ieq]), @view(params.uprimitive[:,:,ieq]), params.visc_coeff[ieq], params.ω, params.mesh, params.basis, params.metrics, params.inputs, iel, ieq, params.QT, SD, params.AD)
        end
        
    end
    
    params.rhs_diff_el .= @views (params.rhs_diffξ_el .+ params.rhs_diffη_el)
end

function _expansion_inviscid!(params, iel, ::CL, QT::Inexact, SD::NSD_1D, AD::FD)

    for ieq = 1:params.neqs
        for i = 1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,1]
            params.RHS[ip,ieq] = 0.0
        end
    end
    
    nothing
end

function _expansion_inviscid!(params, iel, ::CL, QT::Inexact, SD::NSD_1D, AD::ContGal)
    
    for ieq = 1:params.neqs
        for i=1:params.mesh.ngl
            dFdξ = 0.0
            for k = 1:params.mesh.ngl
                dFdξ += params.basis.dψ[k,i]*params.F[k,1,ieq]
            end
            params.rhs_el[iel,i,1,ieq] -= params.ω[i]*dFdξ - params.ω[i]*params.S[i,1,ieq]
        end
    end
end


function _expansion_inviscid!(params, iel, ::CL, QT::Inexact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(params, iel, ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)

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

function _expansion_inviscid!(params, iel, ::CL, QT::Exact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(params, iel, ::CL, QT::Exact, SD::NSD_2D, AD::ContGal)
    
    N = params.mesh.ngl
    Q = N + 1
    for ieq=1:params.neqs
        for l=1:Q
            for k=1:Q
                ωJac = params.ω[k]*params.ω[l]*params.metrics.Je[iel,k,l]
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                for n = 1:N
                    for m = 1:N
                        dFdξ += params.basis.dψ[m,k]* params.basis.ψ[n,l]*params.F[m,n,ieq]
                        dFdη +=  params.basis.ψ[m,k]*params.basis.dψ[n,l]*params.F[m,n,ieq]
                        
                        dGdξ += params.basis.dψ[m,k]* params.basis.ψ[n,l]*params.G[m,n,ieq]
                        dGdη +=  params.basis.ψ[m,k]*params.basis.dψ[n,l]*params.G[m,n,ieq]
                    end
                end
                
                dξdx_kl = params.metrics.dξdx[iel,k,l]
                dξdy_kl = params.metrics.dξdy[iel,k,l]
                dηdx_kl = params.metrics.dηdx[iel,k,l]
                dηdy_kl = params.metrics.dηdy[iel,k,l]
                for j = 1:N
                    for i = 1:N
                        dFdx = dFdξ*dξdx_kl + dFdη*dηdx_kl
                        dGdx = dGdξ*dξdx_kl + dGdη*dηdx_kl

                        dFdy = dFdξ*dξdy_kl + dFdη*dηdy_kl
                        dGdy = dGdξ*dξdy_kl + dGdη*dηdy_kl
                        
                        auxi = ωJac*params.basis.ψ[i,k]*params.basis.ψ[j,l]*((dFdx + dGdy) - params.S[i,j,ieq])
                        params.rhs_el[iel,i,j,ieq] -= auxi
                    end
                end
            end
        end
    end
end

function _expansion_inviscid!(params, iel, ::NCL, QT::Inexact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(params, iel, ::NCL, QT::Inexact, SD::NSD_2D, AD::ContGal)
    
    for ieq=1:params.neqs
        for j=1:params.mesh.ngl
            for i=1:params.mesh.ngl
                ωJac = params.ω[i]*params.ω[j]*params.metrics.Je[iel,i,j]
                
                dFdξ = 0.0; dFdη = 0.0
                dGdξ = 0.0; dGdη = 0.0
                dpdξ = 0.0; dpdη = 0.0               
                for k = 1:params.mesh.ngl
                    dFdξ += params.basis.dψ[k,i]*params.F[k,j,ieq]
                    dFdη += params.basis.dψ[k,j]*params.F[i,k,ieq]
                    
                    dGdξ += params.basis.dψ[k,i]*params.G[k,j,ieq]
                    dGdη += params.basis.dψ[k,j]*params.G[i,k,ieq]
                    
                    dpdξ += params.basis.dψ[k,i]*params.uprimitive[k,j,params.neqs+1]
                    dpdη += params.basis.dψ[k,j]*params.uprimitive[i,k,params.neqs+1]
                end
                dξdx_ij = params.metrics.dξdx[iel,i,j]
                dξdy_ij = params.metrics.dξdy[iel,i,j]
                dηdx_ij = params.metrics.dηdx[iel,i,j]
                dηdy_ij = params.metrics.dηdy[iel,i,j]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij            
                dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij

                dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij            
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
                dpdx = dpdξ*dξdx_ij + dpdη*dηdx_ij            
                dpdy = dpdξ*dξdy_ij + dpdη*dηdy_ij

                ρij = params.uprimitive[i,j,1]
                uij = params.uprimitive[i,j,2]
                vij = params.uprimitive[i,j,3]
                
                if (ieq == 1)
                    auxi = ωJac*(dFdx + dGdy)
                elseif(ieq == 2)
                    auxi = ωJac*(uij*dFdx + vij*dGdy + dpdx/ρij)
                elseif(ieq == 3)
                    auxi = ωJac*(uij*dFdx + vij*dGdy + dpdy/ρij - params.S[i,j,ieq])
                elseif(ieq == 4)
                    auxi = ωJac*(uij*dFdx + vij*dGdy)
                end
                
                params.rhs_el[iel,i,j,ieq] -= auxi
            end
        end
    end        
end


function _expansion_inviscid!(params, iel, ::NCL, QT::Exact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(params, iel, ::NCL, QT::Exact, SD::NSD_2D, AD::ContGal)

    N = params.mesh.ngl
    Q = N + 1

    for l=1:Q
        for k=1:Q
            ωJac = params.ω[k]*params.ω[l]*params.metrics.Je[iel,k,l]
            
            dρudξ = 0.0; dρudη = 0.0
            dρvdξ = 0.0; dρvdη = 0.0
            dudξ = 0.0; dudη = 0.0
            dvdξ = 0.0; dvdη = 0.0
            dθdξ = 0.0; dθdη = 0.0
            dpdξ = 0.0; dpdη = 0.0         
            
            ρkl = 0.0; ukl = 0.0; vkl = 0.0; Skl = 0.0
            for n=1:N
                for m=1:N
                    ψmk = params.basis.ψ[m,k]
                    ψnl = params.basis.ψ[n,l]
                    
                    dψmk_ψnl = params.basis.dψ[m,k]* params.basis.ψ[n,l]
                    ψmk_dψnl = params.basis.ψ[m,k]*params.basis.dψ[n,l]
                    
                    dρudξ += dψmk_ψnl*params.F[m,n,1]
                    dρudη +=  ψmk_dψnl*params.F[m,n,1]
                    
                    dρvdξ += dψmk_ψnl*params.G[m,n,1]
                    dρvdη +=  ψmk_dψnl*params.G[m,n,1]
                    
                    dudξ += dψmk_ψnl*params.uprimitive[m,n,2]
                    dudη +=  ψmk_dψnl*params.uprimitive[m,n,2]

                    dvdξ += dψmk_ψnl*params.uprimitive[m,n,3]
                    dvdη +=  ψmk_dψnl*params.uprimitive[m,n,3]
                    
                    dθdξ += dψmk_ψnl*params.uprimitive[m,n,4]
                    dθdη +=  ψmk_dψnl*params.uprimitive[m,n,4]

                    dpdξ += dψmk_ψnl*params.uprimitive[m,n,params.neqs+1]
                    dpdη +=  ψmk_dψnl*params.uprimitive[m,n,params.neqs+1]

                    ρkl += ψmk*ψnl*params.uprimitive[m,n,1]
                    ukl += ψmk*ψnl*params.uprimitive[m,n,2]
                    vkl += ψmk*ψnl*params.uprimitive[m,n,3]
                    Skl += ψmk*ψnl*params.S[m,n,3]
                end
            end

            dξdx_kl = params.metrics.dξdx[iel,k,l]
            dξdy_kl = params.metrics.dξdy[iel,k,l]
            dηdx_kl = params.metrics.dηdx[iel,k,l]
            dηdy_kl = params.metrics.dηdy[iel,k,l]
            
            dρudx = dρudξ*dξdx_kl + dρudη*dηdx_kl            
            dρudy = dρudξ*dξdy_kl + dρudη*dηdy_kl
            dρvdx = dρvdξ*dξdx_kl + dρvdη*dηdx_kl            
            dρvdy = dρvdξ*dξdy_kl + dρvdη*dηdy_kl
            
            dudx = dudξ*dξdx_kl + dudη*dηdx_kl            
            dudy = dudξ*dξdy_kl + dudη*dηdy_kl
            
            dvdx = dvdξ*dξdx_kl + dvdη*dηdx_kl            
            dvdy = dvdξ*dξdy_kl + dvdη*dηdy_kl
            
            dθdx = dθdξ*dξdx_kl + dθdη*dηdx_kl            
            dθdy = dθdξ*dξdy_kl + dθdη*dηdy_kl

            dpdx = dpdξ*dξdx_kl + dpdη*dηdx_kl            
            dpdy = dpdξ*dξdy_kl + dpdη*dηdy_kl


            for j=1:N
                for i=1:N

                    ψikψjl = params.basis.ψ[i,k]*params.basis.ψ[j,l]
                    
                    params.rhs_el[iel,i,j,1] -= ψikψjl*ωJac*(dρudx + dρvdy)
                    
                    params.rhs_el[iel,i,j,2] -= ψikψjl*ωJac*(ukl*dudx + vkl*dudy + dpdx/ρkl)
                    params.rhs_el[iel,i,j,3] -= ψikψjl*ωJac*(ukl*dvdx + vkl*dvdy + dpdy/ρkl - Skl)
                    params.rhs_el[iel,i,j,4] -= ψikψjl*ωJac*(ukl*dθdx + vkl*dθdy)
                end
            end
            
        end
    end
end


function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeffieq, ω,
                          mesh, basis, metrics, inputs, iel, ieq, QT::Inexact, SD::NSD_2D, ::FD)
    nothing
end

function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeffieq, ω,
                          mesh, basis, metrics, inputs, iel, ieq, QT::Inexact, SD::NSD_2D, ::ContGal)
    
    for l = 1:mesh.ngl
        for k = 1:mesh.ngl
            ωJac = ω[k]*ω[l]*metrics.Je[iel,k,l]
            
            dqdξ = 0.0
            dqdη = 0.0
            @turbo for ii = 1:mesh.ngl
                dqdξ += basis.dψ[ii,k]*uprimitiveieq[ii,l]
                dqdη += basis.dψ[ii,l]*uprimitiveieq[k,ii]
            end
            dξdx_kl = metrics.dξdx[iel,k,l]
            dξdy_kl = metrics.dξdy[iel,k,l]
            dηdx_kl = metrics.dηdx[iel,k,l]
            dηdy_kl = metrics.dηdy[iel,k,l]
            
            auxi = dqdξ*dξdx_kl + dqdη*dηdx_kl
            dqdx = visc_coeffieq*auxi
            
            auxi = dqdξ*dξdy_kl + dqdη*dηdy_kl
            dqdy = visc_coeffieq*auxi
            
            ∇ξ∇u_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
            ∇η∇u_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac     
            
            @turbo for i = 1:mesh.ngl
                dhdξ_ik = basis.dψ[i,k]
                dhdη_il = basis.dψ[i,l]
                
                rhs_diffξ_el[i,l] -= dhdξ_ik * ∇ξ∇u_kl
                rhs_diffη_el[k,i] -= dhdη_il * ∇η∇u_kl
            end
        end
    end  
end

function  _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeff, ω, mesh, basis, metrics, inputs, iel, ieq, QT::Exact, SD::NSD_2D, ::FD)
    nothing
end

function  _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeff, ω, mesh, basis, metrics, inputs, iel, ieq, QT::Exact, SD::NSD_2D, ::ContGal)
    
    N = params.mesh.ngl
    Q = N + 1

    for l=1:Q
        for k=1:Q
            ωJac = params.ω[k]*params.ω[l]*params.metrics.Je[iel,k,l]
            
            dqdξ = 0.0; dqdη = 0.0
            ρkl = 0.0; ukl = 0.0; vkl = 0.0; Skl = 0.0
            for n=1:N
                for m=1:N
                    ψmk = params.basis.ψ[m,k]
                    ψnl = params.basis.ψ[n,l]
                    
                    dψmk_ψnl = params.basis.dψ[m,k]* params.basis.ψ[n,l]
                    ψmk_dψnl = params.basis.ψ[m,k]*params.basis.dψ[n,l]
                    
                    dqdξ += dψmk_ψnl*params.uprimitiveieq[m,n]
                    dqdη += ψmk_dψnl*params.uprimitiveieq[m,n]
                    ukl  +=  ψmk*ψnl*params.uprimitiveieq[m,n]
                    
                end
            end

            dξdx_kl = params.metrics.dξdx[iel,k,l]
            dξdy_kl = params.metrics.dξdy[iel,k,l]
            dηdx_kl = params.metrics.dηdx[iel,k,l]
            dηdy_kl = params.metrics.dηdy[iel,k,l]
            
            dqdx = dqdξ*dξdx_kl + dqdη*dηdx_kl
            dqdx = dqdx*visc_coeff[2]

            dqdy = dqdξ*dξdy_kl + dqdη*dηdy_kl
            dqdy = dqdy*visc_coeff[2]
            
            ∇ξ∇u_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
            ∇η∇u_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac     
            
            ###### W I P ######
            for j=1:N
                for i=1:N

                    dhdξ_ik = basis.dψ[i,k]
                    dhdη_il = basis.dψ[i,l]
                    
                    rhs_diffξ_el[i,l] -= dhdξ_ik * ∇ξ∇u_kl
                    rhs_diffη_el[k,i] -= dhdη_il * ∇η∇u_kl
                    
                    #params.rhs_diffξ_el[iel,i,j,2] -=
                    #params.rhs_diffξ_el[iel,i,j,3] -=
                    #params.rhs_diffξ_el[iel,i,j,4] -=
                    
                    #params.rhs_diffη_el[iel,i,j,2] -=
                    #params.rhs_diffη_el[iel,i,j,3] -=
                    #params.rhs_diffη_el[iel,i,j,4] -=
                end
            end
            
        end
    end
end
