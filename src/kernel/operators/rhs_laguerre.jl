#---------------------------------------------------------------------------
# Optimized (more coud possibly be done)
#---------------------------------------------------------------------------
function build_rhs_laguerre!(RHS, u, params, time)
    #
    # build_rhs()! is called by TimeIntegrators.jl -> time_loop!() via ODEProblem(rhs!, u, tspan, params)
    #
   _build_rhs_laguerre!(RHS, u, params, time)
    
end

function resetRHSToZero_inviscid_laguerre!(params)
    fill!(params.rhs_el_lag, zero(params.T))   
    fill!(params.RHS_lag,    zero(params.T))
end

function resetRHSToZero_viscous_laguerre!(params)
    fill!(params.rhs_diff_el_lag,  zero(params.T))
    fill!(params.rhs_diffξ_el_lag, zero(params.T))
    fill!(params.rhs_diffη_el_lag, zero(params.T))
    fill!(params.RHS_visc_lag,     zero(params.T))
end

function uToPrimitives_laguerre!(neqs, uprimitive, u, uauxe, connijk_lag, ngl, ngr, npoin, δtotal_energy, iel, PT, ::CL, ::AbstractPert, SD::NSD_1D)
   nothing
end


function uToPrimitives_laguerre!(neqs, uprimitive, u, uauxe, connijk_lag, ngl, ngr, npoin, δtotal_energy, iel, PT, ::CL, ::TOTAL, SD::NSD_2D)

    if typeof(PT) == CompEuler

        PhysConst = PhysicalConst{Float64}()

        for j=1:ngr, i=1:ngl

            m1 = connijk_lag[iel,i,j]
            m2 = m1 + npoin
            m3 = m2 + npoin
            m4 = m3 + npoin

            uprimitive[i,j,1] = u[m1]
            uprimitive[i,j,2] = u[m2]/u[m1]
            uprimitive[i,j,3] = u[m3]/u[m1]
            uprimitive[i,j,4] = u[m4]/u[m1] - δtotal_energy*0.5*(uprimitive[i,j,2]^2 + uprimitive[i,j,3]^2)
            #Tracers
            if(neqs > 4)
                mieq = m4
                for ieq = 5:neqs
                    mieq = mieq + npoin
                    uprimitive[i,j,ieq] = u[mieq]
                end
            end
            #Pressure:
            uprimitive[i,j,end] = perfectGasLaw_ρθtoP(PhysConst, ρ=uprimitive[i,j,1], θ=uprimitive[i,j,4])

        end

    elseif typeof(PT) == AdvDiff

        for j=1:ngr, i=1:ngl

            mieq = connijk_lag[iel,i,j]
            uprimitive[i,j,1] = u[mieq]

            for ieq = 2:neqs
                mieq = mieq + npoin
                uprimitive[i,j,ieq] = u[mieq]
            end

        end
    end
end

function uToPrimitives_laguerre!(neqs, uprimitive, u, uauxe, connijk_lag, ngl, ngr, npoin, δtotal_energy, iel, PT, ::CL, ::PERT, SD::NSD_2D)

    if typeof(PT) == CompEuler
        PhysConst = PhysicalConst{Float64}()

        for j=1:ngr, i=1:ngl

            m1 = connijk_lag[iel,i,j]
            m2 = m1 + npoin
            m3 = m2 + npoin
            m4 = m3 + npoin

            uprimitive[i,j,1] = u[m1] + uauxe[m1,1]
            uprimitive[i,j,2] = u[m2]/uprimitive[i,j,1]#(u[m2]+uauxe[m1,2])/uauxe[m1,1]#uprimitive[i,j,1]
            uprimitive[i,j,3] = u[m3]/uprimitive[i,j,1]#(u[m3]+uauxe[m1,3])/uauxe[m1,1]#uprimitive[i,j,1]
            #uprimitive[i,j,4] = (u[m4] + uprimitive[i,j,1]*uauxe[m1,4])/uprimitive[i,j,1] # CHECK THIS FOR ENE- δtotal_energy*0.5*(uprimitive[i,j,2]^2 + uprimitive[i,j,3]^2)
            uprimitive[i,j,4] = (u[m4] + uauxe[m1,4])/uprimitive[i,j,1] - uauxe[m1,4]/uauxe[m1,1]#(u[m4] + uauxe[m1,4])/uprimitive[i,j,1]#(u[m4] + uauxe[m1,4])/uauxe[m1,1] - uauxe[m1,4]/uauxe[m1,1]

            #Tracers
            if(neqs > 4)
                mieq = m4
                for ieq = 5:neqs
                    mieq = mieq + npoin
                    uprimitive[i,j,ieq] = u[mieq] + uauxe[mieq,1]
                end
            end

            #Pressure:
            uprimitive[i,j,end] = perfectGasLaw_ρθtoP(PhysConst, ρ=uprimitive[i,j,1], θ=u[m4]+uauxe[m1,4])
        end

    elseif typeof(PT) == AdvDiff

        for j=1:ngr, i=1:ngl

            mieq = connijk_lag[iel,i,j]
            uprimitive[i,j,1] = u[mieq]

            for ieq = 2:neqs
                mieq = mieq + npoin
                uprimitive[i,j,ieq] = u[mieq]
            end

        end
    end

end


function _build_rhs_laguerre!(RHS, u, params, time)

    T       = Float64
    SD      = params.SD
    QT      = params.QT
    CL      = params.CL
    neqs    = params.neqs
    ngl     = params.mesh.ngl
    nelem   = params.mesh.nelem
    npoin   = params.mesh.npoin
    lsource = params.inputs[:lsource]
    
    #-----------------------------------------------------------------------------------
    # Inviscid rhs:
    #-----------------------------------------------------------------------------------    
    resetRHSToZero_inviscid_laguerre!(params) 
    
    #filter!(u, params, SD)
     
    inviscid_rhs_el_laguerre!(u, params, params.mesh.connijk_lag, params.mesh.x, params.mesh.y, lsource, SD)
    DSS_rhs_laguerre!(params.RHS_lag, params.rhs_el_lag, params.mesh.connijk_lag, params.mesh.nelem_semi_inf, ngl, params.mesh.ngr, neqs, SD, params.AD)
    #@info params.rhs_el_lag[1,1,2,2], params.mesh.connijk_lag[1,1,2]
    #@info params.rhs_el_lag[20,5,2,2], params.mesh.connijk_lag[20,5,2]
    #@info params.RHS_lag[6481,2]
    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)

        resetRHSToZero_viscous_laguerre!(params)
        
        viscous_rhs_el_laguerre!(u, params, SD)
        DSS_rhs_laguerre!(params.RHS_visc_lag, params.rhs_diff_el_lag, params.mesh.connijk_lag, params.mesh.nelem_semi_inf, ngl, params.mesh.ngr, neqs, SD, params.AD)
        
        params.RHS_lag[:,:] .= @view(params.RHS_lag[:,:]) .+ @view(params.RHS_visc_lag[:,:])
    end
    
    for ieq=1:neqs
        divide_by_mass_matrix!(@view(params.RHS_lag[:,ieq]), params.vaux, params.Minv, neqs, npoin, params.AD)
    end
    
end

function inviscid_rhs_el_laguerre!(u, params, connijk_lag, x, y, lsource, SD::NSD_1D)

    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    xmax = params.xmax
    xmin = params.xmin
    ymax = params.ymax   
    for iel=1:params.mesh.nelem_semi_inf

        uToPrimitives_laguerre!(params.neqs, params.uprimitive_lag, u, params.qp.qe, connijk_lag, params.mesh.ngl, params.mesh.ngr, params.mesh.npoin, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)

        for i=1:params.mesh.ngr
            ip = connijk_lag[iel,i,1]
            user_flux!(@view(params.F_lag[i,:]), @view(params.G_lag[i,:]), SD,
                       @view(params.uaux[ip,:]),
                       @view(params.qp.qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs)

            if lsource
                user_source!(@view(params.S_lag[i,:]),
                             @view(params.uaux[ip,:]),
                             @view(params.qp.qe[ip,:]),          #ρref
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE; neqs=params.neqs, x=x[ip],y=y[ip],xmax=xmax,xmin=xmin,ymax=ymax)
            end
        end

        _expansion_inviscid_laguerre!(u, params.neqs, params.mesh.ngr, params.basis_lag.dψ, params.ω_lag, params.F_lag, params.S_lag,
        params.metrics_lag.Je, params.metrics_lag.dξdx, params.rhs_el_lag, iel, params.CL, params.QT, SD, params.AD)

    end
end


function inviscid_rhs_el_laguerre!(u, params, connijk_lag, x, y, lsource, SD::NSD_2D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    xmax = params.xmax
    xmin = params.xmin
    ymax = params.ymax
    for iel=1:params.mesh.nelem_semi_inf
        
        uToPrimitives_laguerre!(params.neqs, params.uprimitive_lag, u, params.qp.qe, connijk_lag, params.mesh.ngl, params.mesh.ngr, params.mesh.npoin, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)       
 
        for j=1:params.mesh.ngr, i=1:params.mesh.ngl
            ip = connijk_lag[iel,i,j]
            #=if (abs(params.mesh.x[ip])>4990) 
              #@info params.mesh.x[ip], params.mesh.y[ip]
              params.uaux[ip,2] = 0.0 
            end=#
            #if (abs(params.mesh.y[ip]) == ymax) params.uaux[ip,3] = 0.0 end 
            user_flux!(@view(params.F_lag[i,j,:]), @view(params.G_lag[i,j,:]), SD,
                        @view(params.uaux[ip,:]), 
                        @view(params.qp.qe[ip,:]),         #pref
                       params.mesh, params.CL, params.SOL_VARS_TYPE; neqs=params.neqs)
            
            if lsource
                user_source!(@view(params.S_lag[i,j,:]),
                             @view(params.uaux[ip,:]),
                             @view(params.qp.qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE; neqs=params.neqs, x=x[ip],y=y[ip],xmax = xmax, xmin=xmin,ymax=ymax)
            end
        end
        _expansion_inviscid_laguerre!(u, params.neqs, params.mesh.ngl, params.mesh.ngr, params.basis.dψ, params.basis_lag.dψ, params.ω, params.ω_lag, params.F_lag, params.G_lag, params.S_lag, 
        params.metrics_lag.Je, params.metrics_lag.dξdx,params.metrics_lag.dξdy, params.metrics_lag.dηdx, params.metrics_lag.dηdy, params.rhs_el_lag, iel, params.CL, params.QT, SD, params.AD)
    end
    uaux2u!(u, @view(params.uaux[:,:]), params.neqs, params.mesh.npoin)
end

function viscous_rhs_el_laguerre!(u, params, SD::NSD_2D)
    
    for iel=1:params.mesh.nelem_semi_inf
        
        uToPrimitives_laguerre!(params.neqs, params.uprimitive_lag, u, params.qp.qe, params.mesh.connijk_lag, params.mesh.ngl, params.mesh.ngr, params.mesh.npoin, params.inputs[:δtotal_energy], iel, params.PT, params.CL, params.SOL_VARS_TYPE, SD)

        for ieq in params.ivisc_equations    
              _expansion_visc_laguerre!(params.rhs_diffξ_el_lag, params.rhs_diffη_el_lag, params.uprimitive_lag, params.visc_coeff, params.ω, params.ω_lag, params.mesh.ngl, params.mesh.ngr, 
            params.basis.ψ, params.basis.dψ, params.basis_lag.ψ, params.basis_lag.dψ, params.metrics_lag.Je, params.metrics_lag.dξdx, params.metrics_lag.dξdy, params.metrics_lag.dηdx, 
            params.metrics_lag.dηdy, params.inputs, iel, ieq, params.QT, SD, params.AD)
        end
        
    end
    
    #@info maximum(params.rhs_diffξ_el_lag[:,:,:,1]),maximum(params.rhs_diffη_el_lag[:,:,:,1])
    params.rhs_diff_el_lag .= @views (params.rhs_diffξ_el_lag .+ params.rhs_diffη_el_lag)
end



function _expansion_inviscid_laguerre!(params, iel, ::CL, QT::Inexact, SD::NSD_1D, ::FD) nothing end

function _expansion_inviscid_laguerre!(u, neqs, ngr, dψ_lag, ω_lag, F_lag, S_lag, Je, dξdx, rhs_el_lag, iel, ::CL, QT::Inexact, SD::NSD_1D, ::ContGal)

    for ieq = 1:neqs
        for i=1:ngr
            dFdξ = 0.0
            for k = 1:ngr
                dFdξ += dψ_lag[k,i]*F_lag[k,ieq]
            end
            rhs_el_lag[iel,i,ieq] -= dFdξ*Je[iel,i]*dξdx[iel,i]*ω_lag[i]  - ω_lag[i]*S_lag[i,ieq]
        end
    end
end

function _expansion_inviscid_laguerre!(params, iel, ::CL, QT::Inexact, SD::NSD_2D, ::FD) nothing end

function _expansion_inviscid_laguerre!(u, neqs, ngl, ngr, dψ, dψ_lag, ω, ω_lag, F_lag, G_lag, S_lag, Je, dξdx, dξdy, dηdx, dηdy, rhs_el_lag, iel, ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)
    for ieq=1:neqs
        for j=1:ngr
            for i=1:ngl
                
                ωJac = ω[i]*ω_lag[j]*Je[iel,i,j]
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                @turbo for k = 1:ngl
                    dFdξ += dψ[k,i]*F_lag[k,j,ieq]
                    
                    dGdξ += dψ[k,i]*G_lag[k,j,ieq]
                end
                @turbo for k = 1:ngr
                    dFdη += dψ_lag[k,j]*F_lag[i,k,ieq]

                    dGdη += dψ_lag[k,j]*G_lag[i,k,ieq]
                end
                dξdx_ij = dξdx[iel,i,j]
                dξdy_ij = dξdy[iel,i,j]
                dηdx_ij = dηdx[iel,i,j]
                dηdy_ij = dηdy[iel,i,j]
                 
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij

                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
		        auxi = ωJac*((dFdx + dGdy) - S_lag[i,j,ieq])
                rhs_el_lag[iel,i,j,ieq] -= auxi
            end
        end
    end
end

function _expansion_inviscid_laguerre!(params, iel, ::CL, QT::Exact, SD::NSD_2D, ::FD) nothing end

function _expansion_inviscid_laguerre!(params, iel, ::CL, QT::Exact, SD::NSD_2D, ::ContGal)
    
    basis1 = params.basis
    basis2 = params.basis_lag
    ω1 = params.ω
    ω2 = params.ω_lag
    N = params.mesh.ngl
    N_lag = params.mesh.ngr
    Q = N + 1
    Q_lag = N_lag + 1
    for ieq=1:params.neqs
        for l=1:Q_lag
            for k=1:Q
                ωJac = ω1[k]*ω2[l]*params.metrics_lag.Je[iel,k,l]
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                for n = 1:N_lag
                    for m = 1:N
                        dFdξ += basis1.dψ[m,k] * basis2.ψ[n,l]*params.F_lag[m,n,ieq]
                        dFdη += basis1.ψ[m,k] * basis2.dψ[n,l]*params.F_lag[m,n,ieq]
                        
                        dGdξ += basis1.dψ[m,k] * basis2.ψ[n,l]*params.G_lag[m,n,ieq]
                        dGdη += basis1.ψ[m,k] * basis2.dψ[n,l]*params.G_lag[m,n,ieq]
                    end
                end
                
                dξdx_kl = params.metrics_lag.dξdx[iel,k,l]
                dξdy_kl = params.metrics_lag.dξdy[iel,k,l]
                dηdx_kl = params.metrics_lag.dηdx[iel,k,l]
                dηdy_kl = params.metrics_lag.dηdy[iel,k,l]
                for j = 1:N_lag
                    for i = 1:N
                        dFdx = dFdξ*dξdx_kl + dFdη*dηdx_kl
                        dGdx = dGdξ*dξdx_kl + dGdη*dηdx_kl

                        dFdy = dFdξ*dξdy_kl + dFdη*dηdy_kl
                        dGdy = dGdξ*dξdy_kl + dGdη*dηdy_kl
                        
                        auxi = ωJac*basis1.ψ[i,k]*basis2.ψ[j,l]*((dFdx + dGdy) - params.S_lag[i,j,ieq])
                        params.rhs_el_lag[iel,i,j,ieq] -= auxi
                    end
                end
            end
        end
    end
end

function _expansion_inviscid_laguerre!(params, iel, ::NCL, QT::Inexact, SD::NSD_2D, ::FD) nothing end

function _expansion_inviscid_laguerre!(params, iel, ::NCL, QT::Inexact, SD::NSD_2D, ::ContGal)
  
    basis1 = params.basis
    basis2 = params.basis_lag
    ω1 = params.ω
    ω2 = params.ω_lag 
    for ieq=1:params.neqs
        for j=1:params.mesh.ngr
            for i=1:params.mesh.ngl
                ωJac = ω1[i]*ω2[j]*params.metrics_lag.Je[iel,i,j]
                
                dFdξ = 0.0; dFdη = 0.0
                dGdξ = 0.0; dGdη = 0.0
                dpdξ = 0.0; dpdη = 0.0               
                for k = 1:params.mesh.ngl
                    dFdξ += basis1.dψ[k,i]*params.F_lag[k,j,ieq]
                    
                    dGdξ += basis1.dψ[k,i]*params.G_lag[k,j,ieq]
                                        
                    dpdξ += basis1.dψ[k,i]*params.uprimitive_lag[k,j,params.neqs+1]
                end
                for k = 1:params.mesh.ngr
                    dFdη += basis2.dψ[k,j]*params.F_lag[i,k,ieq]

                    dGdη += basis2.dψ[k,j]*params.G_lag[i,k,ieq]

                    dpdη += basis2.dψ[k,j]*params.uprimitive_lag[i,k,params.neqs+1]
                end
                dξdx_ij = params.metrics_lag.dξdx[iel,i,j]
                dξdy_ij = params.metrics_lag.dξdy[iel,i,j]
                dηdx_ij = params.metrics_lag.dηdx[iel,i,j]
                dηdy_ij = params.metrics_lag.dηdy[iel,i,j]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij            
                dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij

                dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij            
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
                dpdx = dpdξ*dξdx_ij + dpdη*dηdx_ij            
                dpdy = dpdξ*dξdy_ij + dpdη*dηdy_ij

                ρij = params.uprimitive_lag[i,j,1]
                uij = params.uprimitive_lag[i,j,2]
                vij = params.uprimitive_lag[i,j,3]
                
                if (ieq == 1)
                    auxi = ωJac*(dFdx + dGdy)
                elseif(ieq == 2)
                    auxi = ωJac*(uij*dFdx + vij*dGdy + dpdx/ρij)
                elseif(ieq == 3)
                    auxi = ωJac*(uij*dFdx + vij*dGdy + dpdy/ρij - params.S_lag[i,j,ieq])
                elseif(ieq == 4)
                    auxi = ωJac*(uij*dFdx + vij*dGdy)
                end
                
                params.rhs_el_lag[iel,i,j,ieq] -= auxi
            end
        end
    end        
end


function _expansion_inviscid_laguerre!(params, iel, ::NCL, QT::Exact, SD::NSD_2D, ::FD) nothing end

function _expansion_inviscid_laguerre!(params, iel, ::NCL, QT::Exact, SD::NSD_2D, ::ContGal)

    
    basis1 = params.basis
    basis2 = params.basis_lag
    ω1 = params.ω
    ω2 = params.ω_lag
    N = params.mesh.ngl
    Q = N + 1
    N_lag = params.mesh.ngr
    Q_lag = N_lag + 1
    for l=1:Q_lag
        for k=1:Q
            ωJac = ω1[k]*ω2[l]*params.metrics_lag.Je[iel,k,l]
            
            dρudξ = 0.0; dρudη = 0.0
            dρvdξ = 0.0; dρvdη = 0.0
            dudξ = 0.0; dudη = 0.0
            dvdξ = 0.0; dvdη = 0.0
            dθdξ = 0.0; dθdη = 0.0
            dpdξ = 0.0; dpdη = 0.0
           
            
            ρkl = 0.0; ukl = 0.0; vkl = 0.0; Skl = 0.0
            for n=1:N_lag
                for m=1:N
                    ψmk = basis1.ψ[m,k]
                    ψnl = basis2.ψ[n,l]
                    
                    dψmk_ψnl = basis1.dψ[m,k]*basis2.ψ[n,l]
                    ψmk_dψnl = basis1.ψ[m,k]*basis2.dψ[n,l]
                    
                    dρudξ += dψmk_ψnl*params.F_lag[m,n,1]
                    dρudη +=  ψmk_dψnl*params.F_lag[m,n,1]
                    
                    dρvdξ += dψmk_ψnl*params.G_lag[m,n,1]
                    dρvdη +=  ψmk_dψnl*params.G_lag[m,n,1]
                    
                    dudξ += dψmk_ψnl*params.uprimitive_lag[m,n,2]
                    dudη +=  ψmk_dψnl*params.uprimitive_lag[m,n,2]

                    dvdξ += dψmk_ψnl*params.uprimitive_lag[m,n,3]
                    dvdη +=  ψmk_dψnl*params.uprimitive_lag[m,n,3]
                    
                    dθdξ += dψmk_ψnl*params.uprimitive_lag[m,n,4]
                    dθdη +=  ψmk_dψnl*params.uprimitive_lag[m,n,4]

                    dpdξ += dψmk_ψnl*params.uprimitive_lag[m,n,params.neqs+1]
                    dpdη +=  ψmk_dψnl*params.uprimitive_lag[m,n,params.neqs+1]

                    ρkl += ψmk*ψnl*params.uprimitive_lag[m,n,1]
                    ukl += ψmk*ψnl*params.uprimitive_lag[m,n,2]
                    vkl += ψmk*ψnl*params.uprimitive_lag[m,n,3]
                    Skl += ψmk*ψnl*params.S_lag[m,n,3]
                end
            end

            dξdx_kl = params.metrics_lag.dξdx[iel,k,l]
            dξdy_kl = params.metrics_lag.dξdy[iel,k,l]
            dηdx_kl = params.metrics_lag.dηdx[iel,k,l]
            dηdy_kl = params.metrics_lag.dηdy[iel,k,l]
            
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


            for j=1:N_lag
                for i=1:N

                    ψikψjl = basis1.ψ[i,k]*basis2.ψ[j,l]
                    
                    params.rhs_el_lag[iel,i,j,1] -= ψikψjl*ωJac*(dρudx + dρvdy)
                    
                    params.rhs_el_lag[iel,i,j,2] -= ψikψjl*ωJac*(ukl*dudx + vkl*dudy + dpdx/ρkl)
                    params.rhs_el_lag[iel,i,j,3] -= ψikψjl*ωJac*(ukl*dvdx + vkl*dvdy + dpdy/ρkl - Skl)
                    params.rhs_el_lag[iel,i,j,4] -= ψikψjl*ωJac*(ukl*dθdx + vkl*dθdy)
                end
            end
            
        end
    end
end



function _expansion_visc_laguerre!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq,
                                   visc_coeffieq, ω, ω_lag, mesh, basis, basis_lag,
                                   metrics, metrics_lag, inputs, iel, ieq,
                                   QT::Inexact, SD::NSD_2D, ::FD)
    nothing
end

function _expansion_visc_laguerre!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq,
                                   visc_coeffieq, ω, ω_lag, ngl, ngr, ψ, dψ, ψ_lag, dψ_lag,
                                   Je, dξdx, dξdy, dηdx, dηdy, inputs, iel, ieq,
                                   QT::Inexact, SD::NSD_2D, ::ContGal)
    
    for l = 1:ngr
        for k = 1:ngl
            ωJac = ω[k]*ω_lag[l]*Je[iel,k,l]
            
            dudξ = 0.0
            dudη = 0.0
            @turbo for ii = 1:ngl
                dudξ += dψ[ii,k]*uprimitiveieq[ii,l,ieq]
            end
            @turbo for ii = 1:ngr
                dudη += dψ_lag[ii,l]*uprimitiveieq[k,ii,ieq]
            end
            #@info ieq, dudξ, dudη
            dξdx_kl = dξdx[iel,k,l]
            dξdy_kl = dξdy[iel,k,l]
            dηdx_kl = dηdx[iel,k,l]
            dηdy_kl = dηdy[iel,k,l]
             
            auxi = dudξ*dξdx_kl + dudη*dηdx_kl
            dudx = visc_coeffieq[ieq]*auxi
            
            auxi = dudξ*dξdy_kl + dudη*dηdy_kl
            dudy = visc_coeffieq[ieq]*auxi
            
            ∇ξ∇u_kl = (dξdx_kl*dudx + dξdy_kl*dudy)*ωJac
            ∇η∇u_kl = (dηdx_kl*dudx + dηdy_kl*dudy)*ωJac     
            
            @turbo for i = 1:ngl
                dhdξ_ik = dψ[i,k]
                
                rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * ∇ξ∇u_kl
            end
            @turbo for i = 1:ngr
                dhdη_il = dψ_lag[i,l]

                rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * ∇η∇u_kl
            end
        end
    end  
end

