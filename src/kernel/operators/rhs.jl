#---------------------------------------------------------------------------
# Optimized (more coud possibly be done)
#---------------------------------------------------------------------------
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
    fill!(params.rhs.rhs_el, zero(params.T))   
    fill!(params.rhs.RHS,    zero(params.T))
end

function reset_filters!(params)
    fill!(params.filter.b,      zero(params.T))
    fill!(params.filter.B,      zero(params.T))
end
function reset_laguerre_filters!(params)
    fill!(params.filter.b_lag,      zero(params.T))
    fill!(params.filter.B_lag,      zero(params.T))
end

function resetRHSToZero_viscous!(params, SD::NSD_2D)
    fill!(params.rhs.rhs_diff_el,  zero(params.T))
    fill!(params.rhs.rhs_diffξ_el, zero(params.T))
    fill!(params.rhs.rhs_diffη_el, zero(params.T))
    fill!(params.rhs.RHS_visc,     zero(params.T))
end

function resetRHSToZero_viscous!(params, SD::NSD_3D)
    fill!(params.rhs.rhs_diff_el,  zero(params.T))
    fill!(params.rhs.rhs_diffξ_el, zero(params.T))
    fill!(params.rhs.rhs_diffη_el, zero(params.T))
    fill!(params.rhs.rhs_diffζ_el, zero(params.T))
    fill!(params.rhs.RHS_visc,     zero(params.T))
end


function uToPrimitives!(neqs, uprimitive, u, uauxe, mesh, δtotal_energy, iel, PT, ::CL, ::AbstractPert, SD::NSD_1D)
    nothing
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


function uToPrimitives!(neqs, uprimitive, u, uauxe, mesh, δtotal_energy, iel, PT, ::CL, ::TOTAL, SD::NSD_3D)
    
    if typeof(PT) == CompEuler
        
        PhysConst = PhysicalConst{Float64}()
        
        for k=1:mesh.ngl, j=1:mesh.ngl, i=1:mesh.ngl
            
            m1 = mesh.connijk[iel,i,j,k]
            m2 = m1 + mesh.npoin
            m3 = m2 + mesh.npoin
            m4 = m3 + mesh.npoin
            m5 = m4 + mesh.npoin
            
            uprimitive[i,j,k,1] = u[m1]
            uprimitive[i,j,k,2] = u[m2]/u[m1]
            uprimitive[i,j,k,3] = u[m3]/u[m1]
            uprimitive[i,j,k,4] = u[m4]/u[m1]
            uprimitive[i,j,k,5] = u[m5]/u[m1] - δtotal_energy*0.5*(uprimitive[i,j,k,2]^2 + uprimitive[i,j,k,3]^2 + uprimitive[i,j,k,4]^2)
            #Tracers
            if(neqs > 5)
                mieq = m5
                for ieq = 6:neqs
                    mieq = mieq + mesh.npoin
                    uprimitive[i,j,k,ieq] = u[mieq]
                end
            end
            #Pressure:
            uprimitive[i,j,k,end] = perfectGasLaw_ρθtoP(PhysConst, ρ=uprimitive[i,j,k,1], θ=uprimitive[i,j,k,5])
        end
        
    elseif typeof(PT) == AdvDiff
        
        for k=1:mesh.ngl, j=1:mesh.ngl, i=1:mesh.ngl
            
            mieq = mesh.connijk[iel,i,j,k]
            uprimitive[i,j,k,1] = u[mieq]
            
            for ieq = 2:neqs
                mieq = mieq + mesh.npoin
                uprimitive[i,j,k,ieq] = u[mieq]
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
            uprimitive[i,j,4] = (u[m4] + uauxe[m1,4])/uprimitive[i,j,1] - uauxe[m1,4]/uauxe[m1,1]#(u[m4] + uauxe[m1,4])/uprimitive[i,j,1]

            #Tracers
            if(neqs > 4)
                mieq = m4
                for ieq = 5:neqs
                    mieq = mieq + mesh.npoin
                    uprimitive[i,j,ieq] = u[mieq] + uauxe[mieq,1]
                end
            end
            
            #Pressure:
            #uprimitive[i,j,end] = perfectGasLaw_ρθtoP(PhysConst, ρ=uprimitive[i,j,1], θ=uprimitive[i,j,4])
            
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
    
    _build_rhs!(@view(params.rhs.RHS[:,:]), u, params, time)
    if (params.laguerre) 
        build_rhs_laguerre!(@view(params.rhs_lag.RHS_lag[:,:]), u, params, time)
        params.rhs.RHS .= @views(params.rhs.RHS .+ params.rhs_lag.RHS_lag)
    end
    RHStoDU!(du, @view(params.rhs.RHS[:,:]), params.neqs, params.mesh.npoin)
    
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

    # Filter 
    if (params.inputs[:lfilter])
        reset_filters!(params)
        if (params.laguerre)
            reset_laguerre_filters!(params)
        end
        filter!(u, params, time, SD,params.SOL_VARS_TYPE)
    end
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    apply_boundary_conditions!(u, params.uaux, time, params.qp.qe,
                               params.mesh, params.metrics, params.basis,
                               params.rhs.RHS, params.rhs.rhs_el, params.ubdy,
                               params.ω, neqs, params.inputs, AD, SD)
    
    inviscid_rhs_el!(u, params, lsource, SD)
    DSS_rhs!(@view(params.rhs.RHS[:,:]), (params.rhs.rhs_el), params.mesh, nelem, ngl, neqs, SD, AD)
    
    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)
        
        resetRHSToZero_viscous!(params, SD)
        
        viscous_rhs_el!(u, params, SD)
        
        DSS_rhs!(@view(params.rhs.RHS_visc[:,:]), (params.rhs.rhs_diff_el), params.mesh, nelem, ngl, neqs, SD, AD)
        
        @views params.rhs.RHS[:,:] .= params.rhs.RHS[:,:] .+ params.rhs.RHS_visc[:,:]
    end
    for ieq=1:neqs
        divide_by_mass_matrix!(@view(params.rhs.RHS[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
    end
    
end

function inviscid_rhs_el!(u, params, lsource, SD::NSD_1D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    xmax = params.xmax
    xmin = params.xmin
    ymax = params.ymax    
    for iel=1:params.mesh.nelem

        uToPrimitives!(params.neqs, params.fluxes.uprimitive, u, params.qp.qe, params.mesh, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)
        
        for i=1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,1]
            
            user_flux!(@view(params.fluxes.F[i,:]), @view(params.fluxes.G[i,:]), SD,
                       @view(params.uaux[ip,:]),
                       @view(params.qp.qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip)
            
            if lsource
                user_source!(@view(params.fluxes.S[i,:]),
                             @view(params.uaux[ip,:]),
                             @view(params.qp.qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE; neqs=params.neqs, x=params.mesh.x[ip],y=params.mesh.y[ip],xmax=params.xmax,xmin=params.xmin,ymax=params.ymax)
            end
        end
        
        _expansion_inviscid!(u, params, iel, params.CL, params.QT, SD, params.AD)
        
    end
end

function inviscid_rhs_el!(u, params, lsource, SD::NSD_2D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    for iel = 1:params.mesh.nelem
        
        for j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,j]
            
            user_flux!(@view(params.fluxes.F[i,j,:]), @view(params.fluxes.G[i,j,:]), SD,
                       @view(params.uaux[ip,:]),
                       @view(params.qp.qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip)
            
            if lsource
                user_source!(@view(params.fluxes.S[i,j,:]),
                             @view(params.uaux[ip,:]),
                             @view(params.qp.qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE;
                             neqs=params.neqs,
                             x=params.mesh.x[ip],
                             y=params.mesh.y[ip],
                             xmax=params.xmax,
                             xmin=params.xmin,
                             ymax=params.ymax)
            end
        end

        _expansion_inviscid!(u, params, iel, params.CL, params.QT, SD, params.AD)
        
    end
end


function inviscid_rhs_el!(u, params, lsource, SD::NSD_3D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    for iel = 1:params.mesh.nelem
        
        for k = 1:params.mesh.ngl, j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,j,k]
            
            user_flux!(@view(params.fluxes.F[i,j,k,:]), @view(params.fluxes.G[i,j,k,:]), @view(params.H[i,j,k,:]),
                       @view(params.uaux[ip,:]),
                       @view(params.qp.qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip)
            
            if lsource
                user_source!(@view(params.fluxes.S[i,j,k,:]),
                             @view(params.uaux[ip,:]),
                             @view(params.qp.qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE; neqs=params.neqs)
            end
        end

        _expansion_inviscid!(u, params, iel, params.CL, params.QT, SD, params.AD)
        
    end
end


function viscous_rhs_el!(u, params, SD::NSD_2D)

    for iel=1:params.mesh.nelem        
        uToPrimitives!(params.neqs, params.fluxes.uprimitive, u, params.qp.qe, params.mesh, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)

        _expansion_visc!(u, params, iel, params.CL, params.QT, SD, params.AD)        
    end
    
    params.rhs.rhs_diff_el .= @views (params.rhs.rhs_diffξ_el .+ params.rhs.rhs_diffη_el)
    
end


function viscous_rhs_el!(u, params, SD::NSD_3D)
    
    for iel=1:params.mesh.nelem        
        uToPrimitives!(params.neqs, params.fluxes.uprimitive, u, params.qp.qe, params.mesh, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)

        _expansion_visc!(u, params, iel, params.CL, params.QT, SD, params.AD)        
    end
    
    params.rhs.rhs_diff_el .= @views (params.rhs.rhs_diffξ_el .+ params.rhs.rhs_diffη_el .+ params.rhs.rhs_diffζ_el)
    
end


function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_1D, AD::FD)
    
    for ieq = 1:params.neqs
        for i = 1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,1]
            if (ip < params.mesh.npoin)
                params.rhs.RHS[ip,ieq] = 0.5*(u[ip+1] - u[ip])/(params.mesh.Δx[ip])
            end
        end
    end
    nothing
end


function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_1D, AD::ContGal)
    
    for ieq = 1:params.neqs
        for i=1:params.mesh.ngl
            dFdξ = 0.0
            for k = 1:params.mesh.ngl
                dFdξ += params.basis.dψ[k,i]*params.fluxes.F[k,ieq]
            end
            params.rhs.rhs_el[iel,i,ieq] -= params.ω[i]*dFdξ - params.ω[i]*params.fluxes.S[i,ieq]
        end
    end
end


function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for ieq=1:params.neqs
        for j=1:params.mesh.ngl
            for i=1:params.mesh.ngl
                ωJac = params.ω[i]*params.ω[j]*params.metrics.Je[iel,i,j]
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                @turbo for k = 1:params.mesh.ngl
                    dFdξ += params.basis.dψ[k,i]*params.fluxes.F[k,j,ieq]
                    dFdη += params.basis.dψ[k,j]*params.fluxes.F[i,k,ieq]
                    
                    dGdξ += params.basis.dψ[k,i]*params.fluxes.G[k,j,ieq]
                    dGdη += params.basis.dψ[k,j]*params.fluxes.G[i,k,ieq]
                end
                dξdx_ij = params.metrics.dξdx[iel,i,j]
                dξdy_ij = params.metrics.dξdy[iel,i,j]
                dηdx_ij = params.metrics.dηdx[iel,i,j]
                dηdy_ij = params.metrics.dηdy[iel,i,j]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
                dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij

                dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
                auxi = ωJac*((dFdx + dGdy) - params.fluxes.S[i,j,ieq])
                params.rhs.rhs_el[iel,i,j,ieq] -= auxi
            end
        end
    end
end


function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)

    for ieq=1:params.neqs
        for k=1:params.mesh.ngl
            for j=1:params.mesh.ngl
                for i=1:params.mesh.ngl
                    ωJac = params.ω[i]*params.ω[j]*params.ω[k]*params.metrics.Je[iel,i,j,k]
                    
                    dFdξ = 0.0
                    dFdη = 0.0
                    dFdζ = 0.0
                    
                    dGdξ = 0.0
                    dGdη = 0.0
                    dGdζ = 0.0

                    dHdξ = 0.0
                    dHdη = 0.0
                    dHdζ = 0.0
                    @turbo for m = 1:params.mesh.ngl
                        dFdξ += params.basis.dψ[m,i]*params.fluxes.F[m,j,k,ieq]
                        dFdη += params.basis.dψ[m,j]*params.fluxes.F[i,m,k,ieq]
                        dFdζ += params.basis.dψ[m,k]*params.fluxes.F[i,j,m,ieq]
                        
                        dGdξ += params.basis.dψ[m,i]*params.fluxes.G[m,j,k,ieq]
                        dGdη += params.basis.dψ[m,j]*params.fluxes.G[i,m,k,ieq]
                        dGdζ += params.basis.dψ[m,k]*params.fluxes.G[i,j,m,ieq]
                        
                        dHdξ += params.basis.dψ[m,i]*params.H[m,j,k,ieq]
                        dHdη += params.basis.dψ[m,j]*params.H[i,m,k,ieq]
                        dHdζ += params.basis.dψ[m,k]*params.H[i,j,m,ieq]
                    end
                    dξdx_ij = params.metrics.dξdx[iel,i,j,k]
                    dξdy_ij = params.metrics.dξdy[iel,i,j,k]
                    dξdz_ij = params.metrics.dξdz[iel,i,j,k]
                    
                    dηdx_ij = params.metrics.dηdx[iel,i,j,k]
                    dηdy_ij = params.metrics.dηdy[iel,i,j,k]
                    dηdz_ij = params.metrics.dηdz[iel,i,j,k]

                    dζdx_ij = params.metrics.dζdx[iel,i,j,k]
                    dζdy_ij = params.metrics.dζdy[iel,i,j,k]
                    dζdz_ij = params.metrics.dζdz[iel,i,j,k]
                    
                    dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij + dFdζ*dζdx_ij
                    dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij + dGdζ*dζdx_ij
                    dHdx = dHdξ*dξdx_ij + dHdη*dηdx_ij + dHdζ*dζdx_ij

                    dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij + dFdζ*dζdy_ij
                    dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij + dGdζ*dζdy_ij
                    dHdy = dHdξ*dξdy_ij + dHdη*dηdy_ij + dHdζ*dζdy_ij
                    
                    dFdz = dFdξ*dξdz_ij + dFdη*dηdz_ij + dFdζ*dζdz_ij
                    dGdz = dGdξ*dξdz_ij + dGdη*dηdz_ij + dGdζ*dζdz_ij
                    dHdz = dHdξ*dξdz_ij + dHdη*dηdz_ij + dHdζ*dζdz_ij
                    
                    auxi = ωJac*((dFdx + dGdy + dHdz) - params.fluxes.S[i,j,k,ieq])
                    params.rhs.rhs_el[iel,i,j,k,ieq] -= auxi
                end
            end
        end
    end
end



function _expansion_inviscid!(u, params, iel, ::CL, QT::Exact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(u, params, iel, ::CL, QT::Exact, SD::NSD_2D, AD::ContGal)
    
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
                        dFdξ += params.basis.dψ[m,k]* params.basis.ψ[n,l]*params.fluxes.F[m,n,ieq]
                        dFdη +=  params.basis.ψ[m,k]*params.basis.dψ[n,l]*params.fluxes.F[m,n,ieq]
                        
                        dGdξ += params.basis.dψ[m,k]* params.basis.ψ[n,l]*params.fluxes.G[m,n,ieq]
                        dGdη +=  params.basis.ψ[m,k]*params.basis.dψ[n,l]*params.fluxes.G[m,n,ieq]
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
                        
                        auxi = ωJac*params.basis.ψ[i,k]*params.basis.ψ[j,l]*((dFdx + dGdy) - params.fluxes.S[i,j,ieq])
                        params.rhs.rhs_el[iel,i,j,ieq] -= auxi
                    end
                end
            end
        end
    end
end

function _expansion_inviscid!(u, params, iel, ::NCL, QT::Inexact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(u, params, iel, ::NCL, QT::Inexact, SD::NSD_2D, AD::ContGal)
    
    for ieq=1:params.neqs
        for j=1:params.mesh.ngl
            for i=1:params.mesh.ngl
                ωJac = params.ω[i]*params.ω[j]*params.metrics.Je[iel,i,j]
                
                dFdξ = 0.0; dFdη = 0.0
                dGdξ = 0.0; dGdη = 0.0
                dpdξ = 0.0; dpdη = 0.0               
                for k = 1:params.mesh.ngl
                    dFdξ += params.basis.dψ[k,i]*params.fluxes.F[k,j,ieq]
                    dFdη += params.basis.dψ[k,j]*params.fluxes.F[i,k,ieq]
                    
                    dGdξ += params.basis.dψ[k,i]*params.fluxes.G[k,j,ieq]
                    dGdη += params.basis.dψ[k,j]*params.fluxes.G[i,k,ieq]
                    
                    dpdξ += params.basis.dψ[k,i]*params.fluxes.uprimitive[k,j,params.neqs+1]
                    dpdη += params.basis.dψ[k,j]*params.fluxes.uprimitive[i,k,params.neqs+1]
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

                ρij = params.fluxes.uprimitive[i,j,1]
                uij = params.fluxes.uprimitive[i,j,2]
                vij = params.fluxes.uprimitive[i,j,3]
                
                if (ieq == 1)
                    auxi = ωJac*(dFdx + dGdy)
                elseif(ieq == 2)
                    auxi = ωJac*(uij*dFdx + vij*dGdy + dpdx/ρij)
                elseif(ieq == 3)
                    auxi = ωJac*(uij*dFdx + vij*dGdy + dpdy/ρij - params.fluxes.S[i,j,ieq])
                elseif(ieq == 4)
                    auxi = ωJac*(uij*dFdx + vij*dGdy)
                end
                
                params.rhs.rhs_el[iel,i,j,ieq] -= auxi
            end
        end
    end        
end


function _expansion_inviscid!(u, params, iel, ::NCL, QT::Exact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(u, params, iel, ::NCL, QT::Exact, SD::NSD_2D, AD::ContGal)

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
                    
                    dρudξ += dψmk_ψnl*params.fluxes.F[m,n,1]
                    dρudη +=  ψmk_dψnl*params.fluxes.F[m,n,1]
                    
                    dρvdξ += dψmk_ψnl*params.fluxes.G[m,n,1]
                    dρvdη +=  ψmk_dψnl*params.fluxes.G[m,n,1]
                    
                    dudξ += dψmk_ψnl*params.fluxes.uprimitive[m,n,2]
                    dudη +=  ψmk_dψnl*params.fluxes.uprimitive[m,n,2]

                    dvdξ += dψmk_ψnl*params.fluxes.uprimitive[m,n,3]
                    dvdη +=  ψmk_dψnl*params.fluxes.uprimitive[m,n,3]
                    
                    dθdξ += dψmk_ψnl*params.fluxes.uprimitive[m,n,4]
                    dθdη +=  ψmk_dψnl*params.fluxes.uprimitive[m,n,4]

                    dpdξ += dψmk_ψnl*params.fluxes.uprimitive[m,n,params.neqs+1]
                    dpdη +=  ψmk_dψnl*params.fluxes.uprimitive[m,n,params.neqs+1]

                    ρkl += ψmk*ψnl*params.fluxes.uprimitive[m,n,1]
                    ukl += ψmk*ψnl*params.fluxes.uprimitive[m,n,2]
                    vkl += ψmk*ψnl*params.fluxes.uprimitive[m,n,3]
                    Skl += ψmk*ψnl*params.fluxes.S[m,n,3]
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
                    
                    params.rhs.rhs_el[iel,i,j,1] -= ψikψjl*ωJac*(dρudx + dρvdy)
                    
                    params.rhs.rhs_el[iel,i,j,2] -= ψikψjl*ωJac*(ukl*dudx + vkl*dudy + dpdx/ρkl)
                    params.rhs.rhs_el[iel,i,j,3] -= ψikψjl*ωJac*(ukl*dvdx + vkl*dvdy + dpdy/ρkl - Skl)
                    params.rhs.rhs_el[iel,i,j,4] -= ψikψjl*ωJac*(ukl*dθdx + vkl*dθdy)
                end
            end
            
        end
    end
end


function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeffieq, ω,
                          mesh, basis, metrics, inputs, iel, ieq, QT::Inexact, SD::NSD_2D, ::FD)
    nothing
end

function _expansion_visc!(u, params, iel, ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)
    
    for ieq=1:params.neqs
        for l = 1:params.mesh.ngl
            for k = 1:params.mesh.ngl
                ωJac = params.ω[k]*params.ω[l]*params.metrics.Je[iel,k,l]
                
                dqdξ = 0.0
                dqdη = 0.0
                @turbo for ii = 1:params.mesh.ngl
                    dqdξ += params.basis.dψ[ii,k]*params.fluxes.uprimitive[ii,l,ieq]
                    dqdη += params.basis.dψ[ii,l]*params.fluxes.uprimitive[k,ii,ieq]
                end
                dξdx_kl = params.metrics.dξdx[iel,k,l]
                dξdy_kl = params.metrics.dξdy[iel,k,l]
                dηdx_kl = params.metrics.dηdx[iel,k,l]
                dηdy_kl = params.metrics.dηdy[iel,k,l]
                
                auxi = dqdξ*dξdx_kl + dqdη*dηdx_kl
                dqdx = params.visc_coeff[ieq]*auxi
                
                auxi = dqdξ*dξdy_kl + dqdη*dηdy_kl
                dqdy = params.visc_coeff[ieq]*auxi
                
                ∇ξ∇u_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
                ∇η∇u_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac     
                
                @turbo for i = 1:params.mesh.ngl
                    dhdξ_ik = params.basis.dψ[i,k]
                    dhdη_il = params.basis.dψ[i,l]
                    
                    params.rhs.rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * ∇ξ∇u_kl
                    params.rhs.rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * ∇η∇u_kl
                end
            end
        end  
    end
end


function _expansion_visc!(u, params, iel, ::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)
    
    for ieq=1:params.neqs
        for m = 1:params.mesh.ngl
            for l = 1:params.mesh.ngl
                for k = 1:params.mesh.ngl
                    ωJac = params.ω[k]*params.ω[l]*params.ω[m]*params.metrics.Je[iel,k,l,m]
                    
                    dqdξ = 0.0
                    dqdη = 0.0
                    dqdζ = 0.0
                    @turbo for ii = 1:params.mesh.ngl
                        dqdξ += params.basis.dψ[ii,k]*params.fluxes.uprimitive[ii,l,m,ieq]
                        dqdη += params.basis.dψ[ii,l]*params.fluxes.uprimitive[k,ii,m,ieq]
                        dqdζ += params.basis.dψ[ii,m]*params.fluxes.uprimitive[k,l,ii,ieq]
                    end
                    dξdx_klm = params.metrics.dξdx[iel,k,l,m]
                    dξdy_klm = params.metrics.dξdy[iel,k,l,m]
                    dξdz_klm = params.metrics.dξdz[iel,k,l,m]
                    
                    dηdx_klm = params.metrics.dηdx[iel,k,l,m]
                    dηdy_klm = params.metrics.dηdy[iel,k,l,m]
                    dηdz_klm = params.metrics.dηdz[iel,k,l,m]
                    
                    dζdx_klm = params.metrics.dζdx[iel,k,l,m]
                    dζdy_klm = params.metrics.dζdy[iel,k,l,m]
                    dζdz_klm = params.metrics.dζdz[iel,k,l,m]
                    
                    auxi = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
                    dqdx = params.visc_coeff[ieq]*auxi
                    
                    auxi = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
                    dqdy = params.visc_coeff[ieq]*auxi
                    
                    auxi = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
                    dqdz = params.visc_coeff[ieq]*auxi
                    
                    ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz)*ωJac
                    ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz)*ωJac
                    ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz)*ωJac 
                    
                    @turbo for i = 1:params.mesh.ngl
                        dhdξ_ik = params.basis.dψ[i,k]
                        dhdη_il = params.basis.dψ[i,l]
                        dhdζ_im = params.basis.dψ[i,m]
                        
                        params.rhs.diffξ_el[iel,i,l,m,ieq] -= dhdξ_ik * ∇ξ∇u_klm
                        params.rhs.rhs_diffη_el[iel,k,i,m,ieq] -= dhdη_il * ∇η∇u_klm
                        params.rhs.rhs_diffζ_el[iel,k,l,i,ieq] -= dhdζ_im * ∇ζ∇u_klm
                    end
                end
            end
        end
    end
end

function  _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeff, ω, mesh, basis, metrics, inputs, iel, ieq, QT::Exact, SD::NSD_2D, ::FD)
    nothing
end

function  _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeff, ω, mesh, basis, metrics, inputs, iel, ieq, QT::Exact, SD::NSD_2D, ::ContGal)
    nothing
end
