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


function reset_filters!(params)
    fill!(params.b, zero(params.T))
    fill!(params.B, zero(params.T))
end

function reset_laguerre_filters!(params)
    fill!(params.b_lag, zero(params.T))
    fill!(params.B_lag, zero(params.T))
end

function resetRHSToZero_viscous!(params, SD::NSD_2D)
    fill!(params.rhs_diff_el,  zero(params.T))
    fill!(params.rhs_diffξ_el, zero(params.T))
    fill!(params.rhs_diffη_el, zero(params.T))
    fill!(params.RHS_visc,     zero(params.T))
end

function resetRHSToZero_viscous!(params, SD::NSD_3D)
    fill!(params.rhs_diff_el,  zero(params.T))
    fill!(params.rhs_diffξ_el, zero(params.T))
    fill!(params.rhs_diffη_el, zero(params.T))
    fill!(params.rhs_diffζ_el, zero(params.T))
    fill!(params.RHS_visc,     zero(params.T))
end


function rhs!(du, u, params, time)
    backend = params.inputs[:backend]
    if (backend == CPU())
        build_rhs!(@view(params.RHS[:,:]), u, params, time)
        if (params.laguerre) 
            build_rhs_laguerre!(@view(params.RHS_lag[:,:]), u, params, time)
            params.RHS .= @views(params.RHS .+ params.RHS_lag)
        end
        RHStoDU!(du, @view(params.RHS[:,:]), params.neqs, params.mesh.npoin)
    else
        if (params.SOL_VARS_TYPE == PERT())
            lpert = true
        else
            lpert = false
        end

        if (params.SD == NSD_1D())
            params.RHS .= TFloat(0.0)
            PhysConst = PhysicalConst{TFloat}()

            k1 = utouaux_gpu!(backend)
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.neqs))

            k = _build_rhs_gpu_v0!(backend,(Int64(params.mesh.ngl)))
            k(params.RHS, u, params.uaux, params.qp.qe, params.mesh.x, TFloat(time), params.mesh.connijk , params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu, 
              PhysConst, params.xmax, params.xmin, params.mesh.ngl, params.neqs, lpert, inputs[:lperiodic_1d], params.mesh.npoin_linear, params.mesh.npoin; 
              ndrange = params.mesh.nelem*params.mesh.ngl,workgroupsize = params.mesh.ngl)

            if (params.laguerre)
                params.RHS_lag .= TFloat(0.0)
                k = _build_rhs_gpu_v0!(backend,(Int64(params.mesh.ngr)))
                k(params.RHS, u, params.uaux, params.qp.qe, params.mesh.x, TFloat(time), params.mesh.connijk_lag , params.basis_lag.dψ, params.ω_lag, params.Minv, 
                  params.flux_lag_gpu, params.source_lag_gpu,
                  PhysConst, params.xmax, params.xmin, params.mesh.ngr, params.neqs, lpert, inputs[:lperiodic_1d], params.mesh.npoin_linear, params.mesh.npoin;
                  ndrange = params.mesh.nelem_semi_inf*params.mesh.ngr,workgroupsize = params.mesh.ngr)
                @inbounds  params.RHS .+= params.RHS_lag
            end
            k1 = RHStodu_gpu!(backend)
            k1(params.RHS,du,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl,params.neqs))
        elseif (params.SD == NSD_3D())
            
            params.RHS .= TFloat(0.0)
            PhysConst = PhysicalConst{TFloat}()
            
            k1 = utouaux_gpu!(backend)
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.neqs))
            
            k = apply_boundary_conditions_gpu_3D!(backend)
            k(@view(params.uaux[:,:]), @view(u[:]), params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z, TFloat(time),params.metrics.nx,params.metrics.ny, params.metrics.nz,
              params.mesh.poin_in_bdy_face,params.qbdy_gpu,params.mesh.ngl,TInt(params.neqs), params.mesh.npoin, lpert;
              ndrange = (params.mesh.nfaces_bdy*params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl))
            KernelAbstractions.synchronize(backend)
            
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.neqs))
            
            k = _build_rhs_gpu_3D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
            k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk, params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx, 
              params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz, params.metrics.Je,
              params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu,
              params.mesh.ngl, TInt(params.neqs), PhysConst, params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert;
              ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
            KernelAbstractions.synchronize(backend)
            if (params.inputs[:lvisc])
                params.RHS_visc .= TFloat(0.0)
                params.rhs_diffξ_el .= TFloat(0.0)
                params.rhs_diffη_el .= TFloat(0.0)
                params.rhs_diffζ_el .= TFloat(0.0)
                params.source_gpu .= TFloat(0.0)



                if inputs[:visc_model] == "av" #Default is artificial viscosity with constant coefficient

                    k = _build_rhs_diff_gpu_3D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                    k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el, params.rhs_diffζ_el, params.uaux, params.qp.qe, params.source_gpu, 
                    params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk, 
                    params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, 
                    params.metrics.dζdz, params.metrics.Je, params.basis.dψ, params.ω, params.Minv, params.visc_coeff, params.mesh.ngl, TInt(params.neqs), PhysConst, lpert; 
                    ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))

                elseif inputs[:visc_model] == "smag"
                    k = _build_rhs_diff_gpu_3D_v1!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                    k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el, params.rhs_diffζ_el, params.uaux, params.qp.qe, params.source_gpu,
                    params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk, 
                    params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, 
                    params.metrics.dζdz, params.metrics.Je, params.basis.dψ, params.ω, params.Minv, params.visc_coeff, params.mesh.ngl, TInt(params.neqs), params.mesh.Δeffective_s, PhysConst, lpert; 
                    ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))

                end
                KernelAbstractions.synchronize(backend)
                # KernelAbstractions.synchronize(backend)
                #@info maximum(params.flux_gpu[:,:,:,:,3]), minimum(params.flux_gpu[:,:,:,:,3]), maximum(params.flux_gpu[:,:,:,:,4]), minimum(params.flux_gpu[:,:,:,:,4])
                
                @inbounds params.RHS .+= params.RHS_visc
            end
            KernelAbstractions.synchronize(backend)
        
            k1 = RHStodu_gpu!(backend)
            k1(params.RHS,du,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl,params.neqs))
            
        elseif (params.SD == NSD_2D())
            params.RHS .= TFloat(0.0)
            PhysConst = PhysicalConst{TFloat}()
            k1 = utouaux_gpu!(backend)
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl, params.neqs))
            
            if (params.inputs[:lfilter])
                params.B .= TFloat(0.0)
                kf = filter_gpu_2d!(backend,(Int64(params.mesh.ngl), Int64(params.mesh.ngl)))
                kf(params.uaux, params.qp.qe, params.B, params.fx, params.fy_t, params.metrics.Je, params.ω, params.ω, params.mesh.connijk, params.Minv, 
                   params.mesh.ngl, params.mesh.ngl, params.neqs, lpert;
                   ndrange = (params.mesh.nelem * params.mesh.ngl, params.mesh.ngl), workgroupsize = (params.mesh.ngl, params.mesh.ngl))
                KernelAbstractions.synchronize(backend)
                if (params.laguerre)
                    params.B_lag .= TFloat(0.0)
                    kf = filter_gpu_2d!(backend,(Int64(params.mesh.ngl), Int64(params.mesh.ngr)))
                    kf(params.uaux, params.qp.qe, params.B_lag, params.fx, params.fy_t_lag, params.metrics_lag.Je, 
                       params.ω, params.ω_lag, params.mesh.connijk_lag, params.Minv, params.mesh.ngl, params.mesh.ngr, params.neqs, lpert;
                       ndrange = (params.mesh.nelem_semi_inf * params.mesh.ngl, params.mesh.ngr), workgroupsize = (params.mesh.ngl, params.mesh.ngr))

                    KernelAbstractions.synchronize(backend)

                    params.B .+= params.B_lag
                end
                if (lpert)
                    params.uaux .= params.B
                else
                    params.uaux .= params.B .+ params.qp.qe
                end
                kf = uauxtou_gpu!(backend)
                kf(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl,params.neqs))
                KernelAbstractions.synchronize(backend)
            end
            k = apply_boundary_conditions_gpu!(backend)
            k(@view(params.uaux[:,:]), @view(u[:]), params.qp.qe, params.mesh.x,params.mesh.y,TFloat(time),params.metrics.nx,params.metrics.ny,
              params.mesh.poin_in_bdy_edge,params.qbdy_gpu,params.mesh.ngl,TInt(params.neqs), params.mesh.npoin,lpert;
              ndrange = (params.mesh.nedges_bdy*params.mesh.ngl), workgroupsize = (params.mesh.ngl))
            KernelAbstractions.synchronize(backend)
            if (params.laguerre)

                k = apply_boundary_conditions_lag_gpu!(backend)
                k(@view(params.uaux[:,:]), @view(u[:]), params.qp.qe, params.mesh.x,params.mesh.y,TFloat(time), params.mesh.connijk_lag,
                  params.qbdy_lag_gpu, params.mesh.ngl, params.mesh.ngr, TInt(params.neqs), params.mesh.npoin, params.mesh.nelem_semi_inf, 
                  params.inputs[:lperiodic_laguerre], lpert;
                  ndrange = (params.mesh.nelem_semi_inf*params.mesh.ngl,params.mesh.ngr), workgroupsize = (params.mesh.ngl,params.mesh.ngr))
                KernelAbstractions.synchronize(backend)
            end

            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl,params.neqs))
            k = _build_rhs_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
            k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.connijk, 
              params.metrics.dξdx, params.metrics.dξdy, params.metrics.dηdx, params.metrics.dηdy, params.metrics.Je,
              params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu, params.mesh.ngl, TInt(params.neqs), PhysConst,
              params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, lpert;
              ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl))
            KernelAbstractions.synchronize(backend)
            if (params.laguerre)
                params.RHS_lag .= TFloat(0.0)

                
                k_lag = _build_rhs_lag_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngr)))
                k_lag(params.RHS_lag, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.connijk_lag, params.metrics_lag.dξdx, params.metrics_lag.dξdy,
                      params.metrics_lag.dηdx, params.metrics_lag.dηdy, params.metrics_lag.Je, params.basis.dψ, params.basis_lag.dψ, params.ω,
                      params.ω_lag, params.Minv, params.flux_lag_gpu, params.source_lag_gpu, params.mesh.ngl, params.mesh.ngr, TInt(params.neqs), PhysConst,
                      params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, lpert;
                      ndrange = (params.mesh.nelem_semi_inf*params.mesh.ngl,params.mesh.ngr), workgroupsize = (params.mesh.ngl,params.mesh.ngr))
                KernelAbstractions.synchronize(backend)
                @inbounds params.RHS .+= params.RHS_lag
                if (params.inputs[:lvisc])
                    params.RHS_visc_lag .= TFloat(0.0)
                    params.rhs_diffξ_el_lag .= TFloat(0.0)
                    params.rhs_diffη_el_lag .= TFloat(0.0)
                    params.source_lag_gpu .= TFloat(0.0)

                    k_diff_lag = _build_rhs_visc_lag_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngr)))
                    k_diff_lag(params.RHS_visc_lag, params.rhs_diffξ_el_lag, params.rhs_diffη_el_lag, params.uaux, params.qp.qe, params.source_lag_gpu, params.mesh.x,
                               params.mesh.y, params.mesh.connijk_lag, params.metrics_lag.dξdx, params.metrics_lag.dξdy, params.metrics_lag.dηdx, params.metrics_lag.dηdy,
                               params.metrics_lag.Je, params.basis.dψ, params.basis_lag.dψ, params.ω, params.ω_lag, params.Minv, params.visc_coeff,
                               params.mesh.ngl, params.mesh.ngr, TInt(params.neqs), PhysConst, lpert;
                               ndrange = (params.mesh.nelem_semi_inf*params.mesh.ngl,params.mesh.ngr), workgroupsize = (params.mesh.ngl,params.mesh.ngr))
                    
                    @inbounds params.RHS .+= params.RHS_visc_lag
                    
                end
                
            end

            if (params.inputs[:lvisc])
                params.RHS_visc .= TFloat(0.0)
                params.rhs_diffξ_el .= TFloat(0.0)
                params.rhs_diffη_el .= TFloat(0.0)
                params.source_gpu .= TFloat(0.0)
                
                k = _build_rhs_diff_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el, params.uaux, params.qp.qe, params.source_gpu, params.mesh.x, params.mesh.y, params.mesh.connijk, 
                  params.metrics.dξdx, params.metrics.dξdy, params.metrics.dηdx, params.metrics.dηdy, params.metrics.Je, params.basis.dψ, params.ω, params.Minv, 
                  params.visc_coeff, params.mesh.ngl, TInt(params.neqs), PhysConst, lpert; ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl))
                KernelAbstractions.synchronize(backend)

                @inbounds params.RHS .+= params.RHS_visc
            end
            #@info maximum(params.RHS), maximum(params.RHS_lag), maximum(params.RHS_visc_lag)
            k1 = RHStodu_gpu!(backend)
            k1(params.RHS,du,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl,params.neqs))
            
        end
    end
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
    xmin    = params.mesh.xmin
    xmax    = params.mesh.xmax
    ymin    = params.mesh.ymin
    ymax    = params.mesh.ymax
    zmin    = params.mesh.zmin
    zmax    = params.mesh.zmax    
    
    #-----------------------------------------------------------------------------------
    # Inviscid rhs:
    #-----------------------------------------------------------------------------------    
    resetRHSToZero_inviscid!(params)
    if (params.inputs[:lfilter])
        reset_filters!(params)
        if (params.laguerre)
            reset_laguerre_filters!(params)
            filter!(u, params, time, params.uaux, params.mesh.connijk, params.metrics.Je, SD, params.SOL_VARS_TYPE;
                    connijk_lag = params.mesh.connijk_lag, Je_lag = params.metrics_lag.Je)
        else
            filter!(u, params, time, params.uaux, params.mesh.connijk, params.metrics.Je, SD, params.SOL_VARS_TYPE)
        end
    end
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    apply_boundary_conditions!(u, params.uaux, time, params.qp.qe,
                               params.mesh.x, params.mesh.y, params.mesh.z, params.metrics.nx, params.metrics.ny, params.metrics.nz, params.mesh.npoin, params.mesh.npoin_linear, 
                               params.mesh.poin_in_bdy_edge, params.mesh.poin_in_bdy_face, params.mesh.nedges_bdy, params.mesh.nfaces_bdy, params.mesh.ngl, 
                               params.mesh.ngr, params.mesh.nelem_semi_inf, params.basis.ψ, params.basis.dψ,
                               xmax, ymax, zmax, xmin, ymin, zmin, params.RHS, params.rhs_el, params.ubdy,
                               params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem, params.mesh.bdy_edge_type,
                               params.ω, neqs, params.inputs, AD, SD)
    
    inviscid_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, params.mesh.x, params.mesh.y,lsource, SD)
    
    DSS_rhs!(params.RHS, params.rhs_el, params.mesh.connijk, nelem, ngl, neqs, SD, AD)

    #-----------------------------------------------------------------------------------
    # Kessler:
    #-----------------------------------------------------------------------------------
    #kessler(params.mp, params, @view(params.uaux[:,:]), params.qp.qe, params.mesh)
    
    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)
        
        resetRHSToZero_viscous!(params, SD)
        
        viscous_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, SD)
        
        DSS_rhs!(params.RHS_visc, params.rhs_diff_el, params.mesh.connijk, nelem, ngl, neqs, SD, AD)
        
        params.RHS[:,:] .= @view(params.RHS[:,:]) .+ @view(params.RHS_visc[:,:])
    end
    for ieq=1:neqs
        divide_by_mass_matrix!(@view(params.RHS[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
    end
end

function inviscid_rhs_el!(u, params, connijk, qe, x, y, lsource, SD::NSD_1D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)

    xmin = params.xmin; xmax = params.xmax; ymax = params.ymax
    for iel=1:params.mesh.nelem
        
        for i=1:params.mesh.ngl
            ip = connijk[iel,i,1]
            
            user_primitives!(@view(params.uaux[ip,:]), @view(qe[ip,:]), @view(params.uprimitive[i,:]), params.SOL_VARS_TYPE)

            user_flux!(@view(params.F[i,:]), @view(params.G[i,:]), SD,
                       @view(params.uaux[ip,:]),
                       @view(qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip)
            
            if lsource
                user_source!(@view(params.S[i,:]),
                             @view(params.uaux[ip,:]),
                             @view(qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE;
                             neqs=params.neqs, x=x[ip],y=y[ip],xmax=xmax,xmin=xmin)
            end
        end
        
        _expansion_inviscid!(u, params.neqs, params.mesh.ngl, params.basis.dψ, params.ω, params.F, params.S, params.rhs_el, iel, params.CL, params.QT, SD, params.AD)
        
    end
end

function inviscid_rhs_el!(u, params, connijk, qe, x, y, lsource, SD::NSD_2D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    xmin = params.xmin; xmax = params.xmax; ymax = params.ymax
    for iel = 1:params.mesh.nelem

        for j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = connijk[iel,i,j]
            
            user_primitives!(@view(params.uaux[ip,:]), @view(qe[ip,:]), @view(params.uprimitive[i,j,:]), params.SOL_VARS_TYPE)

            user_flux!(@view(params.F[i,j,:]), @view(params.G[i,j,:]), SD,
                       @view(params.uaux[ip,:]),
                       @view(qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip)
            
            if lsource
                user_source!(@view(params.S[i,j,:]),
                             @view(params.uaux[ip,:]),
                             @view(qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE;
                             neqs=params.neqs, x=x[ip], y=y[ip], xmax=xmax, xmin=xmin, ymax=ymax)
            end
        end

        _expansion_inviscid!(u,
                             params.neqs, params.mesh.ngl,
                             params.basis.dψ, params.ω,
                             params.F, params.G, params.S,
                             params.metrics.Je,
                             params.metrics.dξdx, params.metrics.dξdy,
                             params.metrics.dηdx, params.metrics.dηdy,
                             params.rhs_el, iel, params.CL, params.QT, SD, params.AD)
        
    end
end


function inviscid_rhs_el!(u, params, connijk, qe, x, y, lsource, SD::NSD_3D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    for iel = 1:params.mesh.nelem

        for k = 1:params.mesh.ngl, j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = connijk[iel,i,j,k]
            
            user_flux!(@view(params.F[i,j,k,:]), @view(params.G[i,j,k,:]), @view(params.H[i,j,k,:]),
                       @view(params.uaux[ip,:]),
                       @view(qe[ip,:]),         #pref
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip)
            
            if lsource
                user_source!(@view(params.S[i,j,k,:]),
                             @view(params.uaux[ip,:]),
                             @view(qe[ip,:]),          #ρref 
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE; neqs=params.neqs)
            end
        end

        _expansion_inviscid!(u,
                             params.neqs, params.mesh.ngl,
                             params.basis.dψ, params.ω,
                             params.F, params.G, params.H, params.S,
                             params.metrics.Je,
                             params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz,
                             params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz,
                             params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz,
                             params.rhs_el, iel, params.CL, params.QT, SD, params.AD) 
    end
end


function viscous_rhs_el!(u, params, connijk, qe, SD::NSD_2D)
    
    for iel=1:params.mesh.nelem
        
        for j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = connijk[iel,i,j]

            user_primitives!(@view(params.uaux[ip,:]),@view(qe[ip,:]),@view(params.uprimitive[i,j,:]),params.SOL_VARS_TYPE)
        end

        for ieq in params.ivisc_equations
            _expansion_visc!(params.rhs_diffξ_el,
                             params.rhs_diffη_el,
                             params.uprimitive,
                             params.visc_coeff,
                             params.ω,
                             params.mesh.ngl,
                             params.basis.dψ,
                             params.metrics.Je,
                             params.metrics.dξdx, params.metrics.dξdy,
                             params.metrics.dηdx, params.metrics.dηdy,
                             params.inputs, iel, ieq, params.QT, SD, params.AD)
        end
        
    end
    
    params.rhs_diff_el .= @views (params.rhs_diffξ_el .+ params.rhs_diffη_el)
    
end


function viscous_rhs_el!(u, params, connijk, qe, SD::NSD_3D)
    
    for iel=1:params.mesh.nelem        
        
        for k = 1:params.mesh.ngl, j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = connijk[iel,i,j,k]

            user_primitives!(@view(params.uaux[ip,:]),@view(qe[ip,:]),@view(params.uprimitive[i,j,k,:]),params.SOL_VARS_TYPE)
        end
        for ieq in params.ivisc_equations
            _expansion_visc!(params.rhs_diffξ_el, params.rhs_diffη_el, params.rhs_diffζ_el, params.uprimitive, 
                             params.visc_coeff, params.ω, params.mesh.ngl, params.basis.dψ, params.metrics.Je, params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, 
                             params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx,params.metrics.dζdy, params.metrics.dζdz, params.inputs, iel, ieq, params.QT, SD, params.AD)        
        end
    end
    
    params.rhs_diff_el .= @views (params.rhs_diffξ_el .+ params.rhs_diffη_el .+ params.rhs_diffζ_el)
    
end


function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_1D, AD::FD)
    
    for ieq = 1:params.neqs
        for i = 1:params.mesh.ngl
            ip = params.mesh.connijk[iel,i,1]
            if (ip < params.mesh.npoin)
                params.RHS[ip,ieq] = 0.5*(u[ip+1] - u[ip])/(params.mesh.Δx[ip])
            end
        end
    end
    nothing
end


function _expansion_inviscid!(u, neqs, ngl, dψ, ω, F, S, rhs_el, iel, ::CL, QT::Inexact, SD::NSD_1D, AD::ContGal)
    
    for ieq = 1:neqs
        for i=1:ngl
            dFdξ = 0.0
            for k = 1:ngl
                dFdξ += dψ[k,i]*F[k,ieq]
            end
            rhs_el[iel,i,ieq] -= ω[i]*dFdξ - ω[i]*S[i,ieq]
        end
    end
end


function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(u, neqs, ngl, dψ, ω, F, G, S, Je, dξdx, dξdy, dηdx, dηdy, rhs_el, iel, ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)
    
    for ieq=1:neqs
        for j=1:ngl
            for i=1:ngl
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                @turbo for k = 1:ngl
                    dFdξ += dψ[k,i]*F[k,j,ieq]
                    dFdη += dψ[k,j]*F[i,k,ieq]
                    
                    dGdξ += dψ[k,i]*G[k,j,ieq]
                    dGdη += dψ[k,j]*G[i,k,ieq]
                end
                dξdx_ij = dξdx[iel,i,j]
                dξdy_ij = dξdy[iel,i,j]
                dηdx_ij = dηdx[iel,i,j]
                dηdy_ij = dηdy[iel,i,j]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
                dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij

                dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
                auxi = ωJac*((dFdx + dGdy) - S[i,j,ieq])
                rhs_el[iel,i,j,ieq] -= auxi
            end
        end
    end
end


#function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)
function _expansion_inviscid!(u, neqs, ngl, dψ, ω, F, G, H, S, Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, rhs_el, iel, ::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)
    for ieq=1:neqs
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                    
                    dFdξ = 0.0
                    dFdη = 0.0
                    dFdζ = 0.0
                    
                    dGdξ = 0.0
                    dGdη = 0.0
                    dGdζ = 0.0

                    dHdξ = 0.0
                    dHdη = 0.0
                    dHdζ = 0.0
                    @turbo for m = 1:ngl
                        dFdξ += dψ[m,i]*F[m,j,k,ieq]
                        dFdη += dψ[m,j]*F[i,m,k,ieq]
                        dFdζ += dψ[m,k]*F[i,j,m,ieq]
                        
                        dGdξ += dψ[m,i]*G[m,j,k,ieq]
                        dGdη += dψ[m,j]*G[i,m,k,ieq]
                        dGdζ += dψ[m,k]*G[i,j,m,ieq]
                        
                        dHdξ += dψ[m,i]*H[m,j,k,ieq]
                        dHdη += dψ[m,j]*H[i,m,k,ieq]
                        dHdζ += dψ[m,k]*H[i,j,m,ieq]
                    end
                    dξdx_ij = dξdx[iel,i,j,k]
                    dξdy_ij = dξdy[iel,i,j,k]
                    dξdz_ij = dξdz[iel,i,j,k]
                    
                    dηdx_ij = dηdx[iel,i,j,k]
                    dηdy_ij = dηdy[iel,i,j,k]
                    dηdz_ij = dηdz[iel,i,j,k]

                    dζdx_ij = dζdx[iel,i,j,k]
                    dζdy_ij = dζdy[iel,i,j,k]
                    dζdz_ij = dζdz[iel,i,j,k]
                    
                    dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij + dFdζ*dζdx_ij
                    dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij + dGdζ*dζdx_ij
                    dHdx = dHdξ*dξdx_ij + dHdη*dηdx_ij + dHdζ*dζdx_ij

                    dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij + dFdζ*dζdy_ij
                    dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij + dGdζ*dζdy_ij
                    dHdy = dHdξ*dξdy_ij + dHdη*dηdy_ij + dHdζ*dζdy_ij
                    
                    dFdz = dFdξ*dξdz_ij + dFdη*dηdz_ij + dFdζ*dζdz_ij
                    dGdz = dGdξ*dξdz_ij + dGdη*dηdz_ij + dGdζ*dζdz_ij
                    dHdz = dHdξ*dξdz_ij + dHdη*dηdz_ij + dHdζ*dζdz_ij
                    
                    auxi = ωJac*((dFdx + dGdy + dHdz) - S[i,j,k,ieq])
                    rhs_el[iel,i,j,k,ieq] -= auxi
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
                          ngl, dψ, Je, dξdx, dξdy, dηdx, dηdy, inputs, iel, ieq, QT::Inexact, SD::NSD_2D, ::ContGal)
    
    for l = 1:ngl
        for k = 1:ngl
            ωJac = ω[k]*ω[l]*Je[iel,k,l]
            
            dqdξ = 0.0
            dqdη = 0.0
            @turbo for ii = 1:ngl
                dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,ieq]
                dqdη += dψ[ii,l]*uprimitiveieq[k,ii,ieq]
            end
            dξdx_kl = dξdx[iel,k,l]
            dξdy_kl = dξdy[iel,k,l]
            dηdx_kl = dηdx[iel,k,l]
            dηdy_kl = dηdy[iel,k,l]
            
            auxi = dqdξ*dξdx_kl + dqdη*dηdx_kl
            dqdx = visc_coeffieq[ieq]*auxi
            
            auxi = dqdξ*dξdy_kl + dqdη*dηdy_kl
            dqdy = visc_coeffieq[ieq]*auxi
            
            ∇ξ∇u_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
            ∇η∇u_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac     
            
            @turbo for i = 1:ngl
                dhdξ_ik = dψ[i,k]
                dhdη_il = dψ[i,l]
                
                rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * ∇ξ∇u_kl
                rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * ∇η∇u_kl
            end
        end  
    end
end


function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el, uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, inputs, iel, ieq, QT::Inexact, SD::NSD_3D, ::ContGal)
    
    for m = 1:ngl
        for l = 1:ngl
            for k = 1:ngl
                ωJac = ω[k]*ω[l]*ω[m]*Je[iel,k,l,m]
                
                dqdξ = 0.0
                dqdη = 0.0
                dqdζ = 0.0
                @turbo for ii = 1:ngl
                    dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,m,ieq]
                    dqdη += dψ[ii,l]*uprimitiveieq[k,ii,m,ieq]
                    dqdζ += dψ[ii,m]*uprimitiveieq[k,l,ii,ieq]
                end
                dξdx_klm = dξdx[iel,k,l,m]
                dξdy_klm = dξdy[iel,k,l,m]
                dξdz_klm = dξdz[iel,k,l,m]
                
                dηdx_klm = dηdx[iel,k,l,m]
                dηdy_klm = dηdy[iel,k,l,m]
                dηdz_klm = dηdz[iel,k,l,m]
                
                dζdx_klm = dζdx[iel,k,l,m]
                dζdy_klm = dζdy[iel,k,l,m]
                dζdz_klm = dζdz[iel,k,l,m]
                
                auxi = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
                dqdx = visc_coeffieq[ieq]*auxi
                
                auxi = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
                dqdy = visc_coeffieq[ieq]*auxi
                
                auxi = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
                dqdz = visc_coeffieq[ieq]*auxi
                
                ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz)*ωJac
                ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz)*ωJac
                ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz)*ωJac 
                
                @turbo for i = 1:ngl
                    dhdξ_ik = dψ[i,k]
                    dhdη_il = dψ[i,l]
                    dhdζ_im = dψ[i,m]
                    
                    rhs_diffξ_el[iel,i,l,m,ieq] -= dhdξ_ik * ∇ξ∇u_klm
                    rhs_diffη_el[iel,k,i,m,ieq] -= dhdη_il * ∇η∇u_klm
                    rhs_diffζ_el[iel,k,l,i,ieq] -= dhdζ_im * ∇ζ∇u_klm
                end
            end
        end
    end
end


function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el, uprimitiveieq, visc_coeffieq, ω,
                          mesh, basis, metrics, inputs, iel, ieq, QT::Inexact, SD::NSD_3D, ::ContGal)

    for m = 1:mesh.ngl
        for l = 1:mesh.ngl
            for k = 1:mesh.ngl
                ωJac = ω[k]*ω[l]*ω[m]*metrics.Je[iel,k,l,m]
                
                dqdξ = 0.0
                dqdη = 0.0
                dqdζ = 0.0
                @turbo for ii = 1:mesh.ngl
                    dqdξ += basis.dψ[ii,k]*uprimitiveieq[ii,l,m]
                    dqdη += basis.dψ[ii,l]*uprimitiveieq[k,ii,m]
                    dqdζ += basis.dψ[ii,m]*uprimitiveieq[k,l,ii]
                end
                dξdx_klm = metrics.dξdx[iel,k,l,m]
                dξdy_klm = metrics.dξdy[iel,k,l,m]
                dξdz_klm = metrics.dξdz[iel,k,l,m]
                
                dηdx_klm = metrics.dηdx[iel,k,l,m]
                dηdy_klm = metrics.dηdy[iel,k,l,m]
                dηdz_klm = metrics.dηdz[iel,k,l,m]
                
                dζdx_klm = metrics.dζdx[iel,k,l,m]
                dζdy_klm = metrics.dζdy[iel,k,l,m]
                dζdz_klm = metrics.dζdz[iel,k,l,m]
                
                auxi = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
                dqdx = visc_coeffieq*auxi
                
                auxi = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
                dqdy = visc_coeffieq*auxi
                
                auxi = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
                dqdz = visc_coeffieq*auxi
                
                ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz)*ωJac
                ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz)*ωJac
                ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz)*ωJac 
                
                @turbo for i = 1:mesh.ngl
                    dhdξ_ik = basis.dψ[i,k]
                    dhdη_il = basis.dψ[i,l]
                    dhdζ_im = basis.dψ[i,m]
                    
                    rhs_diffξ_el[i,l,m] -= dhdξ_ik * ∇ξ∇u_klm
                    rhs_diffη_el[k,i,m] -= dhdη_il * ∇η∇u_klm
                    rhs_diffζ_el[k,l,i] -= dhdζ_im * ∇ζ∇u_klm
                end
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
