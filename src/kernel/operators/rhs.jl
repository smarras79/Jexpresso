using Distributions

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

function resetRHSToZero_viscous!(params, SD::NSD_1D)
    fill!(params.rhs_diff_el,  zero(params.T))
    fill!(params.rhs_diffξ_el, zero(params.T))
    fill!(params.RHS_visc,     zero(params.T))
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

function resetbdyfluxToZero!(params)
    fill!(params.F_surf,  zero(params.T))
    fill!(params.S_face,  zero(params.T))
    fill!(params.S_flux,  zero(params.T))
end

function reset∇fToZero!(params, SD::NSD_1D)
    fill!(params.rhs_diff_el,  zero(params.T))
    fill!(params.rhs_diffξ_el, zero(params.T))
    fill!(params.RHS_visc,     zero(params.T))
end

function reset∇fToZero!(params)
    fill!(params.∇f,  zero(params.T))
end

function rhs!(du, u, params, time)
    backend = params.inputs[:backend]
    
    if (backend == CPU())
        _build_rhs!(@view(params.RHS[:,:]), u, params, time)

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
            k(params.RHS, u, params.uaux, params.qp.qe, params.mesh.x, TFloat(time), params.mesh.connijk , params.basis.dψ, params.ω, params.Minv, 
              params.flux_gpu, params.source_gpu, 
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
            MicroConst = MicrophysicalConst{TFloat}()
            k1 = utouaux_gpu!(backend)
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.neqs))
            
            if (params.inputs[:lfilter])
                params.B .= TFloat(0.0)
                kf = filter_gpu_3d!(backend,(Int64(params.mesh.ngl), Int64(params.mesh.ngl), Int64(params.mesh.ngl)))
                kf(@view(params.uaux[:,:]), params.qp.qe, params.B, params.fx, params.fy_t, params.fz_t, params.metrics.Je, params.ω, params.ω, params.ω, params.mesh.connijk, params.Minv,
                   params.mesh.ngl, params.mesh.ngl, params.mesh.ngl, params.neqs, lpert;
                   ndrange = (params.mesh.nelem * params.mesh.ngl, params.mesh.ngl, params.mesh.ngl), workgroupsize = (params.mesh.ngl, params.mesh.ngl, params.mesh.ngl))
                KernelAbstractions.synchronize(backend)
                if (lpert)
                    params.uaux[:,1:params.neqs] .= params.B
                else
                    params.uaux .= params.B .+ params.qp.qe
                end
                kf = uauxtou_gpu!(backend)
                kf(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl,params.neqs))
                KernelAbstractions.synchronize(backend)
            end

            k = apply_boundary_conditions_gpu_3D!(backend)
            k(@view(params.uaux[:,:]), @view(u[:]), params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z, TFloat(time),params.metrics.nx,params.metrics.ny, params.metrics.nz,
              params.mesh.poin_in_bdy_face,params.qbdy_gpu,params.mesh.ngl,TInt(params.neqs), params.mesh.npoin, lpert;
              ndrange = (params.mesh.nfaces_bdy*params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl))
            KernelAbstractions.synchronize(backend)
            
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.neqs))
            
            if (inputs[:lmoist])
                k_moist = do_micro_physics_gpu_3D!(backend)
                k_moist(@view(params.uaux[:,:]), params.qp.qe, params.mp.Tabs, params.mp.qn, params.mp.qi, params.mp.qc,
                        params.mp.qr, params.mp.qs, params.mp.qg, params.mp.Pr, params.mp.Ps, params.mp.Pg,
                        params.mp.S_micro, PhysConst, MicroConst, lpert, params.neqs, params.mesh.npoin, params.mesh.z, params.adjusted, params.Pm; ndrange = (params.mesh.npoin))
                k_precip = _build_precipitation_rhs_gpu_3D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                k_precip(params.RHS, @view(params.uaux[:,:]), params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk,
                         params.metrics.dξdz, params.metrics.dηdz, params.metrics.dζdz, params.metrics.Je,
                         params.basis.dψ, params.ω, params.Minv, params.flux_micro, params.source_micro,
                         params.mesh.ngl, TInt(params.neqs), PhysConst, params.mesh.xmax, params.mesh.xmin,
                         params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert,
                         params.mp.Pr, params.mp.Ps, params.mp.Pg, params.mp.qi, params.mp.qn, params.mp.Tabs, params.mp.S_micro, MicroConst;
                         ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
            end
            KernelAbstractions.synchronize(backend)
            k = _build_rhs_gpu_3D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
            k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk, params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx, 
              params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz, params.metrics.Je,
              params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu,
              params.mesh.ngl, TInt(params.neqs), PhysConst, params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert;
              ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
            if (params.inputs[:case] != "bomex")
                k = _build_rhs_gpu_3D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk, params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx, 
                  params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz, params.metrics.Je,
                  params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu,
                  params.mesh.ngl, TInt(params.neqs), PhysConst, params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert;
                  ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
            else
                k = _build_rhs_gpu_3D_v1!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk, params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx, 
                  params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz, params.metrics.Je,
                  params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu,
                  params.mesh.ngl, TInt(params.neqs), PhysConst, params.thermo_params, params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert;
                  ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
            end

            KernelAbstractions.synchronize(backend)
            if (params.inputs[:lvisc])
                params.RHS_visc     .= TFloat(0.0)
                params.rhs_diffξ_el .= TFloat(0.0)
                params.rhs_diffη_el .= TFloat(0.0)
                params.rhs_diffζ_el .= TFloat(0.0)
                params.source_gpu   .= TFloat(0.0)

                if params.VT == AV() #Default is artificial viscosity with constant coefficient

                    k = _build_rhs_diff_gpu_3D_av!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                    k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el, params.rhs_diffζ_el, params.uaux, params.qp.qe, params.source_gpu, 
                      params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk, 
                      params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, 
                      params.metrics.dζdz, params.metrics.Je, params.basis.dψ, params.ω, params.Minv, params.visc_coeff, params.mesh.ngl, TInt(params.neqs), PhysConst, lpert; 
                      ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))

                elseif params.VT == SMAG()
                    k = _build_rhs_diff_gpu_3D_smag!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                    k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el, params.rhs_diffζ_el, params.uaux, params.qp.qe, params.source_gpu,
                      params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk, 
                      params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, 
                      params.metrics.dζdz, params.metrics.Je, params.basis.dψ, params.ω, params.Minv, params.visc_coeff, params.mesh.ngl, TInt(params.neqs), params.mesh.Δeffective_s, PhysConst, lpert; 
                      ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))

                end
                KernelAbstractions.synchronize(backend)
                if (params.inputs[:case] == "bomex")
                    # param_set = TP.ThermodynamicsParameters(TFloat)
                    k_sa = saturation_adjustment_gpu_3D!(backend)
                    k_sa(params.uaux, params.qp.qe, params.mesh.z, params.mesh.connijk, TInt(params.neqs), params.thermo_params, lpert;
                         ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
                    KernelAbstractions.synchronize(backend)
                    
                    kf = uauxtou_gpu!(backend)
                    kf(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl,params.neqs))
                    KernelAbstractions.synchronize(backend)
                end
                
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
                params.RHS_visc     .= TFloat(0.0)
                params.rhs_diffξ_el .= TFloat(0.0)
                params.rhs_diffη_el .= TFloat(0.0)
                params.source_gpu   .= TFloat(0.0)
                
                k = _build_rhs_diff_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el, params.uaux, params.qp.qe, params.source_gpu, params.mesh.x, params.mesh.y, params.mesh.connijk, 
                  params.metrics.dξdx, params.metrics.dξdy, params.metrics.dηdx, params.metrics.dηdy, params.metrics.Je, params.basis.dψ, params.ω, params.Minv, 
                  params.visc_coeff, params.mesh.ngl, TInt(params.neqs), PhysConst, lpert; ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl))
                KernelAbstractions.synchronize(backend)

                @inbounds params.RHS .+= params.RHS_visc
            end
            #@info maximum(params.RHS), maximum(params.RHS_lag), maximum(params.RHS_visc_lag)
            DSS_global_RHS!(@view(params.RHS[:,:]), params.pM, params.neqs)

            k1 = RHStodu_gpu!(backend)
            k1(params.RHS,du,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.mesh.ngl,params.neqs))

        end
end
end

function _build_rhs!(RHS, u, params, time)

    T       = Float64
    SD      = params.SD
    VT      = params.VT
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

    if SD == NSD_1D()
        comm = MPI.COMM_WORLD
    else
        comm = params.mesh.parts.comm
    end
    mpisize = MPI.Comm_size(comm)
    
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
    
    if inputs[:ladapt] == true
        conformity4ncf_q!(params.uaux, params.pM, SD, QT, params.mesh.connijk, params.mesh, params.Minv, params.metrics.Je, params.ω, AD, neqs, params.interp)
    end
    
    resetbdyfluxToZero!(params)
    apply_boundary_conditions_dirichlet!(u, params.uaux, time, params.qp.qe,
                                         params.mesh.coords,
                                         #params.mesh.x, params.mesh.y, params.mesh.z, 
                                         params.metrics.nx, params.metrics.ny, params.metrics.nz, params.mesh.npoin, params.mesh.npoin_linear, 
                                         params.mesh.poin_in_bdy_edge, params.mesh.poin_in_bdy_face, params.mesh.nedges_bdy, params.mesh.nfaces_bdy, params.mesh.ngl, 
                                         params.mesh.ngr, params.mesh.nelem_semi_inf, params.basis.ψ, params.basis.dψ,
                                         xmax, ymax, zmax, xmin, ymin, zmin, params.RHS, params.rhs_el, params.ubdy,
                                         params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem, 
                                         params.mesh.bdy_edge_type, params.mesh.bdy_face_in_elem, params.mesh.bdy_face_type,
                                         params.mesh.connijk, params.metrics.Jef, params.S_face, 
                                         params.S_flux, params.F_surf, params.M_surf_inv, params.M_edge_inv, params.Minv,
                                         params.mp.Tabs, params.mp.qn,
                                         params.ω, neqs, params.inputs, AD, SD)
    
    if (params.inputs[:lmoist])
        if (SD == NSD_3D())
            do_micro_physics!(params.mp.Tabs, params.mp.qn, params.mp.qc, params.mp.qi, params.mp.qr,
                              params.mp.qs, params.mp.qg, params.mp.Pr, params.mp.Ps, params.mp.Pg, params.mp.S_micro,
                              params.mp.qsatt, params.mesh.npoin, params.uaux, params.mesh.z, params.qp.qe, SD, params.SOL_VARS_TYPE)
        else
            do_micro_physics!(params.mp.Tabs, params.mp.qn, params.mp.qc, params.mp.qi, params.mp.qr,
                              params.mp.qs, params.mp.qg, params.mp.Pr, params.mp.Ps, params.mp.Pg, params.mp.S_micro,
                              params.mp.qsatt, params.mesh.npoin, params.uaux, params.mesh.y, params.qp.qe, SD, params.SOL_VARS_TYPE)
        end
        if (params.inputs[:lprecip])
            compute_precipitation_derivatives!(params.mp.dqpdt, params.mp.dqtdt, params.mp.dhldt, params.mp.Pr, params.mp.Ps,
                                               params.mp.Pg, params.mp.Tabs, params.mp.qi, @view(params.uaux[:,1]), @view(params.qp.qe[:,1]), 
                                               params.mesh.nelem, params.mesh.ngl, params.mesh.connijk, params.H,
                                               params.metrics, params.ω, params.basis.dψ, SD, params.SOL_VARS_TYPE)
            if (SD == NSD_3D())
                params.rhs_el[:,:,:,:,5] .-= params.mp.dhldt
                params.rhs_el[:,:,:,:,6] .+= params.mp.dqtdt
                params.rhs_el[:,:,:,:,7] .+= params.mp.dqpdt
            else
                params.rhs_el[:,:,:,4] .-= params.mp.dhldt
                params.rhs_el[:,:,:,5] .+= params.mp.dqtdt
                params.rhs_el[:,:,:,6] .+= params.mp.dqpdt
            end
        end
        uaux2u!(u, params.uaux, params.neqs, params.mesh.npoin)
    end
    
    if(params.inputs[:lsaturation])
        saturation_adjustment(params.uaux, params.qp.qe, params.mesh.z, params.mesh.connijk, params.mesh.nelem, params.mesh.ngl, neqs, params.thermo_params)
        uaux2u!(u, params.uaux, params.neqs, params.mesh.npoin)
    end
    
    inviscid_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, params.mesh.coords, lsource, SD)
    
    if inputs[:ladapt] == true
        DSS_nc_gather_rhs!(params.RHS, SD, QT, params.rhs_el, params.mesh.connijk, params.mesh.poin_in_edge, params.mesh.non_conforming_facets,
                           params.mesh.non_conforming_facets_parents_ghost, params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost, params.mesh.pgip_owner, ngl-1, neqs, params.interp)
    end
    DSS_rhs!(params.RHS, params.rhs_el, params.mesh.connijk, nelem, ngl, neqs, SD, AD)


    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)
        
        resetRHSToZero_viscous!(params, SD)
        
        viscous_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, SD)
        
        # @info "start DSS_rhs_viscous"
        if inputs[:ladapt] == true
            DSS_nc_gather_rhs!(params.RHS_visc, SD, QT, params.rhs_diff_el, params.mesh.connijk,
                               params.mesh.poin_in_edge, params.mesh.non_conforming_facets,
                               params.mesh.non_conforming_facets_parents_ghost,
                               params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost, params.mesh.pgip_owner,
                               ngl-1, neqs, params.interp)
        end
        DSS_rhs!(params.RHS_visc, params.rhs_diff_el, params.mesh.connijk, nelem, ngl, neqs, SD, AD)
        params.RHS[:,:] .= @view(params.RHS[:,:]) .+ @view(params.RHS_visc[:,:])
    end
    apply_boundary_conditions_neumann!(u, params.uaux, time, params.qp.qe,
                                       params.mesh.coords,
                                       params.metrics.nx, params.metrics.ny, params.metrics.nz,
                                       params.mesh.npoin, params.mesh.npoin_linear,
                                       params.mesh.poin_in_bdy_edge, params.mesh.poin_in_bdy_face, params.mesh.nedges_bdy, params.mesh.nfaces_bdy, params.mesh.ngl,
                                       params.mesh.ngr, params.mesh.nelem_semi_inf, params.basis.ψ, params.basis.dψ,
                                       xmax, ymax, zmax, xmin, ymin, zmin, params.RHS, params.rhs_el, params.ubdy,
                                       params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem, 
                                       params.mesh.bdy_edge_type, params.mesh.bdy_face_in_elem, params.mesh.bdy_face_type,
                                       params.mesh.connijk, params.metrics.Jef, params.S_face, 
                                       params.S_flux, params.F_surf, params.M_surf_inv, params.M_edge_inv, params.Minv,
                                       params.WM.τ_f, params.WM.wθ,
                                       params.mp.Tabs, params.mp.qn,
                                       params.ω, neqs, params.inputs, AD, SD) 
    
    DSS_global_RHS!(@view(params.RHS[:,:]), params.pM, params.neqs)

    for ieq=1:neqs
        divide_by_mass_matrix!(@view(params.RHS[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
        # @info "ieq", ieq
        if inputs[:ladapt] == true
            
            DSS_nc_scatter_rhs!(@view(params.RHS[:,ieq]), SD, QT, selectdim(params.rhs_el, ndims(params.rhs_el), ieq), params.mesh.connijk, params.mesh.poin_in_edge, params.mesh.non_conforming_facets,
                                params.mesh.non_conforming_facets_children_ghost, params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.cgip_ghost, params.mesh.cgip_owner, ngl-1, params.interp)
        end
    end
end

function inviscid_rhs_el!(u, params, connijk, qe, coords, lsource, SD::NSD_1D)
    
    #u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    xmin = params.xmin; xmax = params.xmax; ymax = params.ymax
    for iel=1:params.mesh.nelem
        
        for i=1:params.mesh.ngl
            ip = connijk[iel,i,1]
            
            user_primitives!(@view(params.uaux[ip,:]),
                             @view(qe[ip,:]),
                             @view(params.uprimitive[i,:]),
                             params.SOL_VARS_TYPE)

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
                             neqs=params.neqs, x=coords[ip,1], y=0.0, xmax=xmax,xmin=xmin)
            end
        end
        
        _expansion_inviscid!(u, params.neqs, params.mesh.ngl,
                             params.basis.dψ, params.ω,
                             params.F, params.S,
                             params.rhs_el,
                             iel, params.CL, params.QT, SD, params.AD)
        
    end
end

function inviscid_rhs_el!(u, params, connijk, qe, coords, lsource, SD::NSD_2D)
    
    PhysConst = PhysicalConst{Float64}()

    u_element_wise = zeros(params.mesh.ngl, params.mesh.ngl, params.neqs)

    lkep = false
    
    xmin = params.xmin; xmax = params.xmax; ymax = params.ymax
    for iel = 1:params.mesh.nelem

       
        if lkep
            for j = 1:params.mesh.ngl, i=1:params.mesh.ngl
                ip = connijk[iel,i,j]
                
                user_primitives!(@view(params.uaux[ip,:]),@view(qe[ip,:]),@view(params.uprimitive[i,j,:]), params.SOL_VARS_TYPE)
                
                
                # b. Use the map to find the global point index
                global_point_idx = connijk[iel, i, j]
                
                # c. Find the starting index for this point's data in the flat vector `u`
                start_idx = (global_point_idx - 1) * params.neqs + 1
                
                # d. Copy the 'neqs' variables (ρ, ρu, ρv, E) from u to your 4D array
                u_element_wise[i, j, :] = u[start_idx : start_idx + params.neqs - 1]
            end
        end
        
        for j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = connijk[iel,i,j]
            
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
                             neqs=params.neqs,
                             x=coords[ip,1], y=coords[ip,2],
                             xmax=xmax, xmin=xmin,
                             ymax=ymax)
                
                if (params.inputs[:lmoist])
                    add_micro_precip_sources!(params.mp, params.mp.flux_lw[ip],
                                              params.mp.flux_sw[ip], params.mp.Tabs[ip],
                                              params.mp.S_micro[ip],
                                              @view(params.S[i,j,:]), @view(params.uaux[ip,:]),
                                              params.mp.qn[ip], @view(qe[ip,:]),
                                              SD, params.SOL_VARS_TYPE)
                end
            end

            #=  SM  if luser_function
            user_function!(@view(params.fijk[i,j,:]), SD,
            @view(params.uaux[ip,:]),
            @view(qe[ip,:]),
            params.mesh,
            params.CL, params.SOL_VARS_TYPE;
            neqs=params.neqs, iel=iel, ip=ip)
            end
            =#
        end
        #= SM
        _∇f!(params.∇f_el, params.fijk,
        params.mesh.ngl,
        params.basis.dψ, params.ω,
        params.metrics.Je,
        params.metrics.dξdx, params.metrics.dξdy,
        params.metrics.dηdx, params.metrics.dηdy,
        iel, params.CL, params.QT, SD, params.AD)       
        =#
        
        if lkep
            
            _expansion_inviscid_KEP_twopoint!(u_element_wise,
                                              params.uprimitive,
                                              params.neqs, params.mesh.ngl,
                                              params.basis.dψ, params.ω,
                                              params.F, params.G, params.S,
                                              params.metrics.Je,
                                              params.metrics.dξdx, params.metrics.dξdy,
                                              params.metrics.dηdx, params.metrics.dηdy,
                                              params.rhs_el, iel, params.CL, params.QT, SD, params.AD)
            
        else
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

    #= SM params.rhs_el[:,:,:,2] .-= params.∇f_el[:,:,:,1]
    params.rhs_el[:,:,:,3] .-= params.∇f_el[:,:,:,2]=#

end

function inviscid_rhs_el!(u, params, connijk, qe, coords, lsource, SD::NSD_3D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    xmin = params.xmin; xmax = params.xmax; zmax = params.zmax 
    for iel = 1:params.mesh.nelem

        for k = 1:params.mesh.ngl, j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = connijk[iel,i,j,k]
            
            if !(params.inputs[:lsaturation])
                user_flux!(@view(params.F[i,j,k,:]),
                           @view(params.G[i,j,k,:]),
                           @view(params.H[i,j,k,:]),
                           @view(params.uaux[ip,:]),
                           @view(qe[ip,:]),
                           params.mesh,
                           params.CL, params.SOL_VARS_TYPE;
                           neqs=params.neqs, ip=ip)
            else
                user_flux!(@view(params.F[i,j,k,:]),
                           @view(params.G[i,j,k,:]),
                           @view(params.H[i,j,k,:]),
                           @view(params.uaux[ip,:]),
                           @view(qe[ip,:]),         #pref
                           params.mesh, params.thermo_params,
                           params.CL, params.SOL_VARS_TYPE;
                           neqs=params.neqs, ip=ip,
                           x=coords[ip,1], y=coords[ip,2], z=coords[ip,3])
            end
            
            if lsource
                user_source!(@view(params.S[i,j,k,:]),
                             @view(params.uaux[ip,:]),
                             @view(qe[ip,:]),          #ρref 
                             params.mesh.npoin,
                             params.CL, params.SOL_VARS_TYPE;
                             neqs=params.neqs,
                             x=coords[ip,1], y=coords[ip,2], z=coords[ip,3],
                             xmax=xmax, xmin=xmin, zmax=zmax)
                
                if (params.inputs[:lmoist])
                    add_micro_precip_sources!(params.mp,
                                              params.mp.flux_lw[ip],
                                              params.mp.flux_sw[ip],
                                              params.mp.Tabs[ip],
                                              params.mp.S_micro[ip],
                                              @view(params.S[i,j,k,:]),
                                              @view(params.uaux[ip,:]),
                                              params.mp.qn[ip],
                                              @view(qe[ip,:]),
                                              SD, params.SOL_VARS_TYPE)
                    if (params.inputs[:LST])
                        large_scale_source!(@view(params.uaux[ip,:]),
                                            @view(qe[ip,:]),
                                            @view(params.S[i,j,k,:]), 
                                            params.LST.Rad_cool[ip],
                                            params.LST.T_adv[ip],
                                            params.LST.q_adv[ip],
                                            params.SOL_VARS_TYPE)
                    end
                end

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
                             params.rhs_el, iel, 
                             params.WM.wθ, params.inputs[:lwall_model],
                             params.mesh.connijk,
                             params.mesh.coords,
                             params.mesh.poin_in_bdy_face, params.mesh.elem_to_face, params.mesh.bdy_face_type,
                             params.CL, params.QT, SD, params.AD) 
    end
end



function viscous_rhs_el!(u, params, connijk, qe, SD::NSD_1D)
    
    Δ = params.mesh.Δeffective_l
    
    for iel=1:params.mesh.nelem
        
        for i=1:params.mesh.ngl
            ip = connijk[iel,i]

            user_primitives!(@view(params.uaux[ip,:]), @view(qe[ip,:]), @view(params.uprimitive[i,:]), params.SOL_VARS_TYPE)
        end

        for ieq = 1:params.neqs
            _expansion_visc!(params.rhs_diffξ_el,
                             params.uprimitive,
                             params.visc_coeff,
                             params.ω,
                             params.mesh.ngl,
                             params.basis.dψ,
                             params.metrics.Je,
                             params.metrics.dξdx,
                             params.inputs, params.rhs_el,
                             iel, ieq, params.QT, params.VT, SD, params.AD; Δ=Δ)
        end
        
    end
    
    params.rhs_diff_el .= @views (params.rhs_diffξ_el)
    
end

function viscous_rhs_el!(u, params, connijk, qe, SD::NSD_2D)
    
    Δ = params.mesh.Δeffective_l
      
    for iel=1:params.mesh.nelem
        
        for j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = connijk[iel,i,j]

            user_primitives!(@view(params.uaux[ip,:]),@view(qe[ip,:]),@view(params.uprimitive[i,j,:]), params.SOL_VARS_TYPE)
        end

        # WIP: DSGS
        # compute_viscosity!(params.μsgs, SD, params.PT,
        #                    @view(params.uaux[ip,:]),
        #                    q1, q2,
        #                    rhs, params.Δt, params.mesh, params.metrics)
        
        for ieq = 1:params.neqs
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
                             params.inputs, params.rhs_el,
                             iel, ieq, params.QT, params.VT, SD, params.AD;
                             Δ=Δ)
        end
        
    end
    
    params.rhs_diff_el .= @views (params.rhs_diffξ_el .+ params.rhs_diffη_el)
    
end


function viscous_rhs_el!(u, params, connijk, qe, SD::NSD_3D)
    
    Δ = params.mesh.Δeffective_l
    
    for iel=1:params.mesh.nelem        
        
        for k = 1:params.mesh.ngl, j = 1:params.mesh.ngl, i=1:params.mesh.ngl
            ip = connijk[iel,i,j,k]

            user_primitives!(@view(params.uaux[ip,:]),
                             @view(qe[ip,:]),
                             @view(params.uprimitive[i,j,k,:]),
                             params.SOL_VARS_TYPE)
        end

        
        for ieq = 1:params.neqs
            _expansion_visc!(params.rhs_diffξ_el,
                             params.rhs_diffη_el,
                             params.rhs_diffζ_el,
                             params.uprimitive,
                             params.visc_coeff,
                             params.ω,
                             params.mesh.ngl,
                             params.basis.dψ,
                             params.metrics.Je,
                             params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, 
                             params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz,
                             params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz,
                             params.inputs, params.rhs_el, iel, ieq,
                             params.WM.τ_f, params.WM.wθ, params.inputs[:lwall_model], params.mesh.connijk,
                             params.mesh.coords,                             
                             params.mesh.poin_in_bdy_face, params.mesh.elem_to_face, params.mesh.bdy_face_type,
                             params.QT, params.VT, SD, params.AD; Δ=Δ)
            
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


function _expansion_inviscid!(u, neqs, ngl,
                              dψ, ω,
                              F, S,
                              rhs_el,
                              iel, ::CL, QT::Inexact, SD::NSD_1D, AD::ContGal)
    
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

using StaticArrays

function _expansion_inviscid_KEP!(u, uprimitive,
                                  neqs, ngl, dψ, ω,
                                  F, G, S,
                                  Je,
                                  dξdx, dξdy,
                                  dηdx, dηdy,
                                  rhs_el, iel,
                                  ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)

    # Temporary array to store the divergence at each quadrature point
    # Using StaticArrays for performance (avoids heap allocation in the loop)
    Div = MVector{4, Float64}(0.0, 0.0, 0.0, 0.0)
    
    # Loop over quadrature points in the element
    for j = 1:ngl
        for i = 1:ngl
            ωJac = ω[i] * ω[j] * Je[iel, i, j]
            
            # --- 1. Compute standard divergence for all equations at point (i,j) ---
            # This is necessary because the KEP form for momentum/energy
            # depends on the divergence of the continuity equation.
            for ieq = 1:neqs
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                @turbo for k = 1:ngl
                    dFdξ += dψ[k, i] * F[k, j, ieq]
                    dFdη += dψ[k, j] * F[i, k, ieq]
                    
                    dGdξ += dψ[k, i] * G[k, j, ieq]
                    dGdη += dψ[k, j] * G[i, k, ieq]
                end
                
                dξdx_ij = dξdx[iel, i, j]
                dξdy_ij = dξdy[iel, i, j]
                dηdx_ij = dηdx[iel, i, j]
                dηdy_ij = dηdy[iel, i, j]

                dFdx = dFdξ * dξdx_ij + dFdη * dηdx_ij
                dGdy = dGdξ * dξdy_ij + dGdη * dηdy_ij
                
                #dFdy = dFdξ * dξdy_ij + dFdη * dηdy_ij                
                #dGdx = dGdξ * dξdx_ij + dGdη * dηdx_ij
                
                Div[ieq] = dFdx + dGdy
            end

            # --- 2. Get primitive variables at the quadrature point (i,j) ---
            # This assumes that the solution points and quadrature points are the same
            # (i.e., a collocation-based method on Gauss-Lobatto nodes).
            ρ  = uprimitive[i, j, 1]
            ρu = uprimitive[i, j, 2]
            ρv = uprimitive[i, j, 3]
            
            inv_ρ = 1.0 / ρ
            u_vel = ρu * inv_ρ
            v_vel = ρv * inv_ρ
            
            # --- 3. Apply the KEP split-form correction ---
            # The split form is: 0.5 * [ (Standard Divergence) + (Mass Flux Divergence) * (Velocity) ]
            div_mass_flux = Div[1]
            
            # x-momentum (ieq=2)
            div_mom_x = 0.5 * (Div[2] + div_mass_flux * u_vel)
            
            # y-momentum (ieq=3)
            div_mom_y = 0.5 * (Div[3] + div_mass_flux * v_vel)
            
            # Energy (ieq=4)
            # This is one possible KEP extension to the energy equation.
            # Other, more complex forms exist, often coupled with entropy stability.
            KE = 0.5 * (u_vel^2 + v_vel^2)
            div_energy = 0.5 * (Div[4] + div_mass_flux * KE)
            
            # --- 4. Update the element-local RHS with the KEP-corrected values ---
            rhs_el[iel, i, j, 1] -= ωJac * (div_mass_flux - S[i, j, 1])
            rhs_el[iel, i, j, 2] -= ωJac * (div_mom_x     - S[i, j, 2])
            rhs_el[iel, i, j, 3] -= ωJac * (div_mom_y     - S[i, j, 3])
            rhs_el[iel, i, j, 4] -= ωJac * (div_energy    - S[i, j, 4])
        end
    end
end

function _expansion_inviscid!(u, neqs, ngl, dψ, ω,
                              F, G, S,
                              Je,
                              dξdx, dξdy,
                              dηdx, dηdy,
                              rhs_el, iel,
                              ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)
    
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
                #dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij

                #dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
                rhs_el[iel,i,j,ieq] -=  ωJac*((dFdx + dGdy) - S[i,j,ieq])
            end
        end
    end
end

function _expansion_inviscid!(u, neqs, ngl, dψ, ω,
                              F, G, H, S,
                              Je,
                              dξdx, dξdy, dξdz,
                              dηdx, dηdy, dηdz,
                              dζdx, dζdy, dζdz,
                              rhs_el, iel,
                              wθ, lwall_model,
                              connijk,
                              coords,
                              poin_in_bdy_face, elem_to_face, bdy_face_type,
                              ::CL, QT::Inexact, SD::NSD_3D, AD::ContGal)
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


function _expansion_visc!(rhs_diffξ_el, uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je, dξdx, inputs, rhs_el, iel, ieq,
                          QT::Inexact, VT::AV, SD::NSD_1D, ::ContGal; Δ=1.0)

    for k = 1:ngl
        ωJac = ω[k]*Je[iel,k]
        
        dqdξ = 0.0
        @turbo for ii = 1:ngl
            dqdξ += dψ[ii,k]*uprimitiveieq[ii,ieq]
        end

        dξdx_kl = dqdξ*dξdx[iel,k]
        dqdx = visc_coeffieq[ieq]*dξdx_kl
        
        ∇ξ∇u_kl = dξdx_kl*dqdx*ωJac
        
        @turbo for i = 1:ngl
            dhdξ_ik = dψ[i,k]
            
            rhs_diffξ_el[iel,i,ieq] -= dhdξ_ik * ∇ξ∇u_kl
        end
    end
end


function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeffieq, ω,
                          mesh, basis, metrics, inputs, rhs_el, iel, ieq,
                          QT::Inexact, VT, SD::NSD_2D, ::FD)
    nothing
end

function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, ieq,
                          QT::Inexact, VT::AV, SD::NSD_2D, ::ContGal; Δ=1.0)
    
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

#
# RHS with SMAG 2D
#
function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, ieq,
                          QT::Inexact, VT::SMAG, SD::NSD_2D, ::ContGal; Δ=1.0)
    
    #
    # Constants for Richardson stability correction
    #
    PhysConst  = PhysicalConst{Float32}()
    Pr_t       = PhysConst.Pr_t
    #
    # Neutral/unstable: Pr_t ≈ 0.7 - 0.85
    # Stable:           Pr_t ≈ 1.0 - 2.0 (usually handled with Richardson corrections)
    # Very unstable:    Pr_t ≈ 1/3
    #
    κ          = PhysConst.κ
    cp         = PhysConst.cp
    C_s        = PhysConst.C_s
    C_s2       = C_s^2
    
    for l = 1:ngl
        for k = 1:ngl
            ωJac = ω[k]*ω[l]*Je[iel,k,l]

            # Quantities for Smagorinsky 
            dudξ = 0.0; dudη = 0.0
            dvdξ = 0.0; dvdη = 0.0
            @turbo for ii = 1:ngl
                dudξ += dψ[ii,k]*uprimitiveieq[ii,l,2]
                dudη += dψ[ii,l]*uprimitiveieq[k,ii,2]

                dvdξ += dψ[ii,k]*uprimitiveieq[ii,l,3]
                dvdη += dψ[ii,l]*uprimitiveieq[k,ii,3]
            end
            dξdx_kl = dξdx[iel,k,l]
            dξdy_kl = dξdy[iel,k,l]
            dηdx_kl = dηdx[iel,k,l]
            dηdy_kl = dηdy[iel,k,l]

            #u
            dudx = dudξ*dξdx_kl + dudη*dηdx_kl
            dudy = dudξ*dξdy_kl + dudη*dηdy_kl
            
            #v
            dvdx = dvdξ*dξdx_kl + dvdη*dηdx_kl
            dvdy = dvdξ*dξdy_kl + dvdη*dηdy_kl
            
            # Smagorinsky
            # Strain rate tensor (symmetric part of velocity gradient)
            S11 = dudx
            S22 = dvdy
            S12 = 0.5 * (dudy + dvdx)
            
            # CORRECT Strain rate magnitude calculation
            # |S| = sqrt(2 * S_ij * S_ij)
            S_ij_S_ij = S11*S11 + S22*S22 + 2.0 * S12*S12
            Sij = sqrt(2.0 * S_ij_S_ij)
            
            # Filter width calculation
            Δ2      = Δ * Δ
            
            # Base Smagorinsky eddy viscosity
            ν_t_base = C_s2 * Δ2 * Sij
            ν_t = ν_t_base
            
            # END Smagorinsky

            # Compute scalar gradient for diffusion iequation by iequation
            dqdξ = 0.0; dqdη = 0.0
            @turbo for ii = 1:ngl
                dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,ieq]
                dqdη += dψ[ii,l]*uprimitiveieq[k,ii,ieq]
            end
            # Transform scalar gradient to physical coordinates
            dqdx_phys = dqdξ*dξdx_kl + dqdη*dηdx_kl
            dqdy_phys = dqdξ*dξdy_kl + dqdη*dηdy_kl

            
            # Determine effective diffusivity based on scalar type
            # TODO: Replace this logic with proper equation identification
            # Common orderings:
            # - Conservative: [ρ, ρu, ρv, ρw, ρE] or [ρ, ρu, ρv, ρw, ρE, ρθ]
            # - Primitive: [ρ, u, v, w, T] or [ρ, u, v, w, p, θ]
            if ieq == 4  # Assuming potential temperature equation is at index 5
                # For temperature: use thermal diffusivity (ν_t / Pr_t)
                ρ           = uprimitiveieq[k,l,1]                
                α_molecular = κ / (ρ * cp)  # Molecular thermal diffusivity
                α_turbulent = ν_t / Pr_t    # Turbulent thermal diffusivity
                effective_diffusivity = visc_coeffieq[ieq] * α_turbulent
                
            else
                # For momentum equations: use momentum diffusivity
                effective_diffusivity = visc_coeffieq[ieq] * ν_t
            end
                        
            # Apply effective diffusivity to scalar gradients
            dqdx = effective_diffusivity * dqdx_phys
            dqdy = effective_diffusivity * dqdy_phys
            
            ∇ξ∇q_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
            ∇η∇q_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac     
            
            @turbo for i = 1:ngl
                dhdξ_ik = dψ[i,k]
                dhdη_il = dψ[i,l]
                
                rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * ∇ξ∇q_kl
                rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * ∇η∇q_kl
            end
        end  
    end
end

function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, ieq,
                          QT::Inexact, VT::SMAG, SD::NSD_2D, ::ContGal; Δ=1.0)
    
    PhysConst  = PhysicalConst{Float32}()
    Pr_t       = PhysConst.Pr_t
    κ          = PhysConst.κ
    cp         = PhysConst.cp
    
    # Test filter width (typically 2× grid filter)
    Δ_test = 2.0 * Δ
    α = Δ_test / Δ  # Filter ratio
    
    eps_dynamic = 1.0e-14  # Prevent division by zero
    
    for l = 1:ngl
        for k = 1:ngl
            ωJac = ω[k]*ω[l]*Je[iel,k,l]

            # ========================================
            # GRID-SCALE VELOCITY GRADIENTS
            # ========================================
            dudξ = 0.0; dudη = 0.0
            dvdξ = 0.0; dvdη = 0.0
            @turbo for ii = 1:ngl
                dudξ += dψ[ii,k]*uprimitiveieq[ii,l,2]
                dudη += dψ[ii,l]*uprimitiveieq[k,ii,2]
                dvdξ += dψ[ii,k]*uprimitiveieq[ii,l,3]
                dvdη += dψ[ii,l]*uprimitiveieq[k,ii,3]
            end
            
            dξdx_kl = dξdx[iel,k,l]
            dξdy_kl = dξdy[iel,k,l]
            dηdx_kl = dηdx[iel,k,l]
            dηdy_kl = dηdy[iel,k,l]

            # Transform to physical space
            dudx = dudξ*dξdx_kl + dudη*dηdx_kl
            dudy = dudξ*dξdy_kl + dudη*dηdy_kl
            dvdx = dvdξ*dξdx_kl + dvdη*dηdx_kl
            dvdy = dvdξ*dξdy_kl + dvdη*dηdy_kl
            
            # Grid-scale strain rate tensor
            S11 = dudx
            S22 = dvdy
            S12 = 0.5 * (dudy + dvdx)
            
            # Grid-scale strain rate magnitude
            S_ij_S_ij = S11*S11 + S22*S22 + 2.0 * S12*S12
            Sij = sqrt(2.0 * S_ij_S_ij)
            
            # ========================================
            # TEST-FILTERED VELOCITIES
            # ========================================
            # Apply simple box filter to velocities at test scale
            # For spectral elements, approximate with local averaging
            u_kl = uprimitiveieq[k,l,2]
            v_kl = uprimitiveieq[k,l,3]
            
            # Test-filtered velocities (simple box filter approximation)
            u_test = u_kl
            v_test = v_kl
            n_avg = 0
            for jj = max(1,k-1):min(ngl,k+1)
                for ll = max(1,l-1):min(ngl,l+1)
                    u_test += uprimitiveieq[jj,ll,2]
                    v_test += uprimitiveieq[jj,ll,3]
                    n_avg += 1
                end
            end
            u_test /= (n_avg + 1)
            v_test /= (n_avg + 1)
            
            # Test-filtered velocity gradients
            # Simplified: use grid gradients for test filter (approximation)
            dudx_test = dudx
            dudy_test = dudy
            dvdx_test = dvdx
            dvdy_test = dvdy
            
            # Test-scale strain rate tensor
            S11_test = dudx_test
            S22_test = dvdy_test
            S12_test = 0.5 * (dudy_test + dvdx_test)
            
            S_test_ij_S_test_ij = S11_test*S11_test + S22_test*S22_test + 2.0*S12_test*S12_test
            Sij_test = sqrt(2.0 * S_test_ij_S_test_ij)
            
            # ========================================
            # GERMANO IDENTITY TERMS
            # ========================================
            # Leonard stresses: L_ij = test_filter(u_i*u_j) - test_filter(u_i)*test_filter(u_j)
            # Approximate with resolved velocities
            L11 = u_kl*u_kl - u_test*u_test
            L12 = u_kl*v_kl - u_test*v_test
            L22 = v_kl*v_kl - v_test*v_test
            
            # M_ij terms for dynamic procedure
            # M_ij = -2*(Δ_test² |S_test| S_test_ij - α² Δ² |S| S_ij)
            Δ2 = Δ * Δ
            Δ_test2 = Δ_test * Δ_test
            
            M11 = -2.0 * (Δ_test2 * Sij_test * S11_test - α*α * Δ2 * Sij * S11)
            M12 = -2.0 * (Δ_test2 * Sij_test * S12_test - α*α * Δ2 * Sij * S12)
            M22 = -2.0 * (Δ_test2 * Sij_test * S22_test - α*α * Δ2 * Sij * S22)
            
            # Compute C_s² using least-squares (Lilly's contraction)
            # C_s² = <L_ij M_ij> / <M_ij M_ij>
            LM = L11*M11 + 2.0*L12*M12 + L22*M22
            MM = M11*M11 + 2.0*M12*M12 + M22*M22
            
            # Dynamic coefficient
            if MM > eps_dynamic
                C_s2_dynamic = LM / MM
            else
                C_s2_dynamic = 0.0
            end
            
            # Clip negative values (optional: allow backscatter by removing this)
            C_s2_dynamic = max(0.0, C_s2_dynamic)
            
            # Upper bound for stability (optional)
            C_s2_max = 0.09  # (C_s_max ≈ 0.3)²
            C_s2_dynamic = min(C_s2_dynamic, C_s2_max)
            
            # ========================================
            # DYNAMIC SMAGORINSKY EDDY VISCOSITY
            # ========================================
            ν_t = C_s2_dynamic * Δ2 * Sij
            
            # ========================================
            # SCALAR DIFFUSION
            # ========================================
            dqdξ = 0.0; dqdη = 0.0
            @turbo for ii = 1:ngl
                dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,ieq]
                dqdη += dψ[ii,l]*uprimitiveieq[k,ii,ieq]
            end
            
            dqdx_phys = dqdξ*dξdx_kl + dqdη*dηdx_kl
            dqdy_phys = dqdξ*dξdy_kl + dqdη*dηdy_kl

            # Effective diffusivity
            if ieq == 4  # Potential temperature
                α_turbulent = ν_t / Pr_t
                effective_diffusivity = visc_coeffieq[ieq] * α_turbulent
            else  # Momentum
                effective_diffusivity = visc_coeffieq[ieq] * ν_t
            end
                        
            dqdx = effective_diffusivity * dqdx_phys
            dqdy = effective_diffusivity * dqdy_phys
            
            ∇ξ∇q_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
            ∇η∇q_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac     
            
            @turbo for i = 1:ngl
                dhdξ_ik = dψ[i,k]
                dhdη_il = dψ[i,l]
                
                rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * ∇ξ∇q_kl
                rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * ∇η∇q_kl
            end
        end  
    end
end

#
# RHS with Wall-Adapting Local Eddy-viscosity (WALE) 2D
#
function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, ieq,
                          QT::Inexact, VT::WALE, SD::NSD_2D, ::ContGal; Δ=1.0)

    #
    # Constants:
    #
    PhysConst = PhysicalConst{Float32}()
    Pr_t      = PhysConst.Pr_t
    κ         = PhysConst.κ
    cp        = PhysConst.cp
    C_w       = 0.5  # WALE constant, typically ~0.5
    C_w2      = C_w^2
    
    for l = 1:ngl
        for k = 1:ngl
            ωJac = ω[k]*ω[l]*Je[iel,k,l]

            # Quantities for Smagorinsky 
            dudξ = 0.0; dudη = 0.0
            dvdξ = 0.0; dvdη = 0.0
            @turbo for ii = 1:ngl
                dudξ += dψ[ii,k]*uprimitiveieq[ii,l,2]
                dudη += dψ[ii,l]*uprimitiveieq[k,ii,2]

                dvdξ += dψ[ii,k]*uprimitiveieq[ii,l,3]
                dvdη += dψ[ii,l]*uprimitiveieq[k,ii,3]
            end
            dξdx_kl = dξdx[iel,k,l]
            dξdy_kl = dξdy[iel,k,l]
            dηdx_kl = dηdx[iel,k,l]
            dηdy_kl = dηdy[iel,k,l]

            #u
            dudx = dudξ*dξdx_kl + dudη*dηdx_kl
            dudy = dudξ*dξdy_kl + dudη*dηdy_kl
            
            #v
            dvdx = dvdξ*dξdx_kl + dvdη*dηdx_kl
            dvdy = dvdξ*dξdy_kl + dvdη*dηdy_kl

            
            #------------------------------------------------------------------------
            # WALE Model
            #------------------------------------------------------------------------
            # 1. Strain-rate tensor (S_ij)
            S11 = dudx
            S22 = dvdy
            S12 = 0.5 * (dudy + dvdx)
            
            # 2. Rotation-rate tensor (Omega_ij)
            O12 = 0.5 * (dudy - dvdx)

            # 3. Square of the velocity gradient tensor (g_sq_ij = g_ik * g_kj)
            #    For 2D: g = [[dudx, dudy], [dvdx, dvdy]]
            g_sq_11 = dudx * dudx + dudy * dvdx
            g_sq_12 = dudx * dudy + dudy * dvdy
            g_sq_21 = dvdx * dudx + dvdy * dvdx
            g_sq_22 = dvdx * dudy + dvdy * dvdy

            # 4. Traceless symmetric part of g_sq (S_d_ij)
            #    Symmetric part of g_sq
            S_g_sq_12 = 0.5 * (g_sq_12 + g_sq_21)
            #    Trace of g_sq (for 2D, S33 is zero)
            trace_g_sq = g_sq_11 + g_sq_22
            #    Compute S_d_ij for 2D
            Sd11 = g_sq_11 - (1.0 / 2.0) * trace_g_sq # Note: 1/2 for 2D, 1/3 for 3D
            Sd22 = g_sq_22 - (1.0 / 2.0) * trace_g_sq
            Sd12 = S_g_sq_12
            
            # 5. Calculate scalar invariants
            S_ij_S_ij   = S11^2 + S22^2 + 2.0 * S12^2
            Sd_ij_Sd_ij = Sd11^2 + Sd22^2 + 2.0 * Sd12^2
            
            # 6. Filter width calculation
            Δ2      = Δ * Δ
            
            # 7. Calculate WALE eddy viscosity (ν_t)
            epsilon = 1.0e-10 # To prevent division by zero
            
            numerator   = (Sd_ij_Sd_ij)^1.5
            denominator = (S_ij_S_ij)^2.5 + (Sd_ij_Sd_ij)^1.25
            
            ν_t = C_w2 * Δ2 * (numerator / (denominator + epsilon))
            #------------------------------------------------------------------------
            # END WALE Model
            #------------------------------------------------------------------------
       
            # Compute scalar gradient for diffusion iequation by iequation
            dqdξ = 0.0; dqdη = 0.0
            @turbo for ii = 1:ngl
                dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,ieq]
                dqdη += dψ[ii,l]*uprimitiveieq[k,ii,ieq]
            end
            # Transform scalar gradient to physical coordinates
            dqdx_phys = dqdξ*dξdx_kl + dqdη*dηdx_kl
            dqdy_phys = dqdξ*dξdy_kl + dqdη*dηdy_kl

            
            # Determine effective diffusivity based on scalar type
            # TODO: Replace this logic with proper equation identification
            # Common orderings:
            # - Conservative: [ρ, ρu, ρv, ρw, ρE] or [ρ, ρu, ρv, ρw, ρE, ρθ]
            # - Primitive: [ρ, u, v, w, T] or [ρ, u, v, w, p, θ]
             # Determine effective diffusivity based on scalar type
            if ieq == 4 # Assuming potential temperature equation
                ρ           = uprimitiveieq[k, l, 1]
                α_turbulent = ν_t / Pr_t
                effective_diffusivity = visc_coeffieq[ieq] * α_turbulent
            else # For momentum equations
                effective_diffusivity = visc_coeffieq[ieq] * ν_t
            end
            
            # Apply effective diffusivity to scalar gradients
            dqdx = effective_diffusivity * dqdx_phys
            dqdy = effective_diffusivity * dqdy_phys
            
            ∇ξ∇q_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
            ∇η∇q_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac     
            
            @turbo for i = 1:ngl
                dhdξ_ik = dψ[i,k]
                dhdη_il = dψ[i,l]
                
                rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * ∇ξ∇q_kl
                rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * ∇η∇q_kl
            end
        end  
    end
end

#
# RHS with Vreman 2D
#
function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, ieq,
                          QT::Inexact, VT::VREM, SD::NSD_2D, ::ContGal; Δ=1.0)

    PhysConst  = PhysicalConst{Float32}()
    Pr_t       = PhysConst.Pr_t
    #
    # Neutral/unstable: Pr_t ≈ 0.7 - 0.85
    # Stable:           Pr_t ≈ 1.0 - 2.0 (usually handled with Richardson corrections)
    # Very unstable:    Pr_t ≈ 1/3
    #
    κ          = PhysConst.κ
    cp         = PhysConst.cp
    C_s        = PhysConst.C_s
    C_s2       = C_s^2
    C_vrem     = 2.5 * C_s2  # Vreman coefficient
    eps_vreman = 1.0e-14  # Safety epsilon
    
    for l = 1:ngl
        for k = 1:ngl
            Je_kl = Je[iel,k,l]
            ωJac = ω[k]*ω[l]*Je_kl

            # Filter width calculation (isotropic)
            Δ2      = Δ * Δ
            
            # Velocity gradients in computational space
            dudξ = 0.0; dudη = 0.0
            dvdξ = 0.0; dvdη = 0.0
            @turbo for ii = 1:ngl
                dudξ += dψ[ii,k]*uprimitiveieq[ii,l,2]
                dudη += dψ[ii,l]*uprimitiveieq[k,ii,2]
                dvdξ += dψ[ii,k]*uprimitiveieq[ii,l,3]
                dvdη += dψ[ii,l]*uprimitiveieq[k,ii,3]
            end
            
            # Metric terms
            dξdx_kl = dξdx[iel,k,l]
            dξdy_kl = dξdy[iel,k,l]
            dηdx_kl = dηdx[iel,k,l]
            dηdy_kl = dηdy[iel,k,l]

            # Transform to physical space
            u11 = dudξ*dξdx_kl + dudη*dηdx_kl  # dudx
            u12 = dudξ*dξdy_kl + dudη*dηdy_kl  # dudy
            u21 = dvdξ*dξdx_kl + dvdη*dηdx_kl  # dvdx
            u22 = dvdξ*dξdy_kl + dvdη*dηdy_kl  # dvdy

            # Vreman β tensor
            β11 = Δ2 * (u11*u11 + u12*u12)
            β12 = Δ2 * (u11*u21 + u12*u22)
            β22 = Δ2 * (u21*u21 + u22*u22)

            B_β = β11*β22 - β12*β12
            
            # Frobenius norm squared of velocity gradient
            u_ij_u_ij = u11*u11 + u12*u12 + u21*u21 + u22*u22
            
            # Vreman eddy viscosity with safety checks
            # At the top, after line 13
           
            if u_ij_u_ij > eps_vreman && B_β > 0.0
                ν_t = C_vrem * sqrt(B_β / u_ij_u_ij)
            else
                ν_t = 0.0
            end
            
            # Scalar gradient
            dqdξ = 0.0; dqdη = 0.0
            @turbo for ii = 1:ngl
                dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,ieq]
                dqdη += dψ[ii,l]*uprimitiveieq[k,ii,ieq]
            end
            
            dqdx_phys = dqdξ*dξdx_kl + dqdη*dηdx_kl
            dqdy_phys = dqdξ*dξdy_kl + dqdη*dηdy_kl
            
            # Effective diffusivity
            if ieq == 4
                ρ           = uprimitiveieq[k,l,1]                
                #α_molecular = κ / (ρ * cp)  # Molecular thermal diffusivity
                α_turbulent = ν_t / Pr_t    # Turbulent thermal diffusivity

                effective_diffusivity = visc_coeffieq[ieq] * α_turbulent
            else
                effective_diffusivity = visc_coeffieq[ieq] * ν_t
            end

            #=if iel == 1 && k == 1 && l == 1 && ieq == 2
                @show C_s
                @show Pr_t
                @show effective_diffusivity
                @show α_turbulent
                @show Δ
                @show ν_t
                @show sqrt(B_β / u_ij_u_ij)  # Before multiplying by C_vrem
            end=#
            
            dqdx = effective_diffusivity * dqdx_phys
            dqdy = effective_diffusivity * dqdy_phys
            
            ∇ξ∇q_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
            ∇η∇q_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac     
            
            @turbo for i = 1:ngl
                dhdξ_ik = dψ[i,k]
                dhdη_il = dψ[i,l]
                
                rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * ∇ξ∇q_kl
                rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * ∇η∇q_kl
            end
        end  
    end
end
#
# END RHS with Vreman 2D
#

function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy, dξdz,
                          dηdx, dηdy, dηdz,
                          dζdx, dζdy, dζdz,
                          inputs,
                          rhs_el,
                          iel, ieq,
                          τ_f, wθ, lwall_model,
                          connijk,
                          coords, 
                          poin_in_bdy_face, elem_to_face, bdy_face_type,
                          QT::Inexact, VT::AV, SD::NSD_3D, ::ContGal; Δ=1.0)

    PhysConst = PhysicalConst{Float32}()
    MPConst   = MicrophysicalConst{Float32}()
        
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



function  _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy, dξdz,
                          dηdx, dηdy, dηdz,
                          dζdx, dζdy, dζdz,
                          inputs,
                          rhs_el,
                          iel, ieq,
                          τ_f, wθ, lwall_model,
                          connijk,
                          coords, 
                          poin_in_bdy_face, elem_to_face, bdy_face_type,
                          QT::Inexact, VT::VREM, SD::NSD_3D, ::ContGal; Δ=1.0)
    
    PhysConst  = PhysicalConst{Float32}()
    Pr_t       = PhysConst.Pr_t
    κ          = PhysConst.κ
    cp         = PhysConst.cp
    C_s        = PhysConst.C_s
    C_s2       = C_s^2
    C_vrem     = 2.5 * C_s2  # Vreman coefficient
    eps_vreman = 1.0e-14  # Safety epsilon
    
    for m = 1:ngl
        for l = 1:ngl
            for k = 1:ngl
                Je_klm = Je[iel,k,l,m]
                ωJac = ω[k]*ω[l]*ω[m]*Je_klm
                
                # Filter width calculation (isotropic)
                Δ2      = Δ * Δ
            #    @info Δ
                # Velocity gradients in computational space
                dudξ = 0.0; dudη = 0.0; dudζ = 0.0
                dvdξ = 0.0; dvdη = 0.0; dvdζ = 0.0
                dwdξ = 0.0; dwdη = 0.0; dwdζ = 0.0
                @turbo for ii = 1:ngl
                    dudξ += dψ[ii,k]*uprimitiveieq[ii,l,m,2]
                    dudη += dψ[ii,l]*uprimitiveieq[k,ii,m,2]
                    dudζ += dψ[ii,m]*uprimitiveieq[k,l,ii,2]
                    
                    dvdξ += dψ[ii,k]*uprimitiveieq[ii,l,m,3]
                    dvdη += dψ[ii,l]*uprimitiveieq[k,ii,m,3]
                    dvdζ += dψ[ii,m]*uprimitiveieq[k,l,ii,3]
                    
                    dwdξ += dψ[ii,k]*uprimitiveieq[ii,l,m,4]
                    dwdη += dψ[ii,l]*uprimitiveieq[k,ii,m,4]
                    dwdζ += dψ[ii,m]*uprimitiveieq[k,l,ii,4]
                end
                
                # Metric terms
                dξdx_klm = dξdx[iel,k,l,m]
                dξdy_klm = dξdy[iel,k,l,m]
                dξdz_klm = dξdz[iel,k,l,m]
                dηdx_klm = dηdx[iel,k,l,m]
                dηdy_klm = dηdy[iel,k,l,m]
                dηdz_klm = dηdz[iel,k,l,m]
                dζdx_klm = dζdx[iel,k,l,m]
                dζdy_klm = dζdy[iel,k,l,m]
                dζdz_klm = dζdz[iel,k,l,m]
                
                # Transform to physical space - 3x3 velocity gradient tensor
                u11 = dudξ*dξdx_klm + dudη*dηdx_klm + dudζ*dζdx_klm  # ∂u/∂x
                u12 = dudξ*dξdy_klm + dudη*dηdy_klm + dudζ*dζdy_klm  # ∂u/∂y
                u13 = dudξ*dξdz_klm + dudη*dηdz_klm + dudζ*dζdz_klm  # ∂u/∂z
                
                u21 = dvdξ*dξdx_klm + dvdη*dηdx_klm + dvdζ*dζdx_klm  # ∂v/∂x
                u22 = dvdξ*dξdy_klm + dvdη*dηdy_klm + dvdζ*dζdy_klm  # ∂v/∂y
                u23 = dvdξ*dξdz_klm + dvdη*dηdz_klm + dvdζ*dζdz_klm  # ∂v/∂z
                
                u31 = dwdξ*dξdx_klm + dwdη*dηdx_klm + dwdζ*dζdx_klm  # ∂w/∂x
                u32 = dwdξ*dξdy_klm + dwdη*dηdy_klm + dwdζ*dζdy_klm  # ∂w/∂y
                u33 = dwdξ*dξdz_klm + dwdη*dηdz_klm + dwdζ*dζdz_klm  # ∂w/∂z
                
                # Vreman β tensor (3D)
                # β_ij = Δ_m^2 * u_im * u_jm (sum over m=1,2,3)
                β11 = Δ2 * (u11*u11 + u12*u12 + u13*u13)
                β12 = Δ2 * (u11*u21 + u12*u22 + u13*u23)
                β13 = Δ2 * (u11*u31 + u12*u32 + u13*u33)
                β22 = Δ2 * (u21*u21 + u22*u22 + u23*u23)
                β23 = Δ2 * (u21*u31 + u22*u32 + u23*u33)
                β33 = Δ2 * (u31*u31 + u32*u32 + u33*u33)
                
                # B_β for 3D
                B_β = β11*β22 + β11*β33 + β22*β33 - (β12*β12 + β13*β13 + β23*β23)
                
                # Frobenius norm squared of 3x3 velocity gradient tensor
                u_ij_u_ij = u11*u11 + u12*u12 + u13*u13 +
                            u21*u21 + u22*u22 + u23*u23 +
                            u31*u31 + u32*u32 + u33*u33
                
                # Vreman eddy viscosity with safety checks
                if u_ij_u_ij > eps_vreman && B_β > 0.0
                    ν_t = C_vrem * sqrt(B_β / u_ij_u_ij)
                else
                    ν_t = 0.0
                end
                
                # Scalar gradient in computational space
                dqdξ = 0.0; dqdη = 0.0; dqdζ = 0.0
                @turbo for ii = 1:ngl
                    dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,m,ieq]
                    dqdη += dψ[ii,l]*uprimitiveieq[k,ii,m,ieq]
                    dqdζ += dψ[ii,m]*uprimitiveieq[k,l,ii,ieq]
                end
                
                # Transform scalar gradient to physical space
                dqdx_phys = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
                dqdy_phys = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
                dqdz_phys = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
                
                # Effective diffusivity
                if ieq == 5  # Energy equation (shifted by 1 due to w at index 4)
                    ρ           = uprimitiveieq[k,l,m,1]                
                    α_molecular = κ / (ρ * cp)  # Molecular thermal diffusivity
                    α_turbulent = ν_t / Pr_t    # Turbulent thermal diffusivity
                    
                    effective_diffusivity = visc_coeffieq[ieq] * α_turbulent
                    #effective_diffusivity = visc_coeffieq[ieq] * ρ * cp * (α_molecular + α_turbulent)
                else
                    effective_diffusivity = visc_coeffieq[ieq] * ν_t
                end
                
                # Apply diffusivity to scalar gradients
                dqdx = effective_diffusivity * dqdx_phys
                dqdy = effective_diffusivity * dqdy_phys
                dqdz = effective_diffusivity * dqdz_phys
                
                # Weak form contributions
                ∇ξ∇q_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz)*ωJac
                ∇η∇q_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz)*ωJac
                ∇ζ∇q_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz)*ωJac
                
                @turbo for i = 1:ngl
                    dhdξ_ik = dψ[i,k]
                    dhdη_il = dψ[i,l]
                    dhdζ_im = dψ[i,m]
                    
                    rhs_diffξ_el[iel,i,l,m,ieq] -= dhdξ_ik * ∇ξ∇q_klm
                    rhs_diffη_el[iel,k,i,m,ieq] -= dhdη_il * ∇η∇q_klm
                    rhs_diffζ_el[iel,k,l,i,ieq] -= dhdζ_im * ∇ζ∇q_klm
                end
            end
        end  
    end
end
               
function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                          uprimitive, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy, dξdz,
                          dηdx, dηdy, dηdz,
                          dζdx, dζdy, dζdz,
                          inputs,
                          rhs_el,
                          iel, ieq,
                          τ_f, wθ, lwall_model,
                          connijk,
                          coords,
                          poin_in_bdy_face, elem_to_face, bdy_face_type, 
                          QT::Inexact, VT::SMAG, SD::NSD_3D, ::ContGal; Δ=1.0)

    #
    # Constants for Richardson stability correction
    #
    PhysConst = PhysicalConst{Float32}()
    Pr_t      = PhysConst.Pr_t            # Turbulent Prandtl number
    g         = PhysConst.g               # Gravitational acceleration (m/s²)
    Ri_crit   = PhysConst.Ri_crit         # Critical Richardson number
    C_s       = PhysConst.C_s             # Smagorinsky constant (typical range: 0.1-0.2)
    
    for m = 1:ngl
        for l = 1:ngl
            for k = 1:ngl
                ωJac = ω[k]*ω[l]*ω[m]*Je[iel,k,l,m]
                
                # Initialize velocity gradients
                dudξ = 0.0; dudη = 0.0; dudζ = 0.0
                dvdξ = 0.0; dvdη = 0.0; dvdζ = 0.0
                dwdξ = 0.0; dwdη = 0.0; dwdζ = 0.0
                
                # Initialize potential temperature gradients
                dθdξ = 0.0; dθdη = 0.0; dθdζ = 0.0

                @turbo for ii = 1:ngl
                    # Velocity gradients
                    dudξ += dψ[ii,k]*uprimitive[ii,l,m,2]
                    dudη += dψ[ii,l]*uprimitive[k,ii,m,2]
                    dudζ += dψ[ii,m]*uprimitive[k,l,ii,2]

                    dvdξ += dψ[ii,k]*uprimitive[ii,l,m,3]
                    dvdη += dψ[ii,l]*uprimitive[k,ii,m,3]
                    dvdζ += dψ[ii,m]*uprimitive[k,l,ii,3]

                    dwdξ += dψ[ii,k]*uprimitive[ii,l,m,4]
                    dwdη += dψ[ii,l]*uprimitive[k,ii,m,4]
                    dwdζ += dψ[ii,m]*uprimitive[k,l,ii,4]
                    
                    # Potential temperature gradients
                    dθdξ += dψ[ii,k]*uprimitive[ii,l,m,5]
                    dθdη += dψ[ii,l]*uprimitive[k,ii,m,5]
                    dθdζ += dψ[ii,m]*uprimitive[k,l,ii,5]
                end
                
                # Cache metric derivatives for efficiency
                dξdx_klm = dξdx[iel,k,l,m]
                dξdy_klm = dξdy[iel,k,l,m]
                dξdz_klm = dξdz[iel,k,l,m]
                
                dηdx_klm = dηdx[iel,k,l,m]
                dηdy_klm = dηdy[iel,k,l,m]
                dηdz_klm = dηdz[iel,k,l,m]
                
                dζdx_klm = dζdx[iel,k,l,m]
                dζdy_klm = dζdy[iel,k,l,m]
                dζdz_klm = dζdz[iel,k,l,m]

                # Transform velocity derivatives to physical coordinates
                dudx = dudξ*dξdx_klm + dudη*dηdx_klm + dudζ*dζdx_klm
                dvdx = dvdξ*dξdx_klm + dvdη*dηdx_klm + dvdζ*dζdx_klm
                dwdx = dwdξ*dξdx_klm + dwdη*dηdx_klm + dwdζ*dζdx_klm
                
                dudy = dudξ*dξdy_klm + dudη*dηdy_klm + dudζ*dζdy_klm
                dvdy = dvdξ*dξdy_klm + dvdη*dηdy_klm + dvdζ*dζdy_klm
                dwdy = dwdξ*dξdy_klm + dwdη*dηdy_klm + dwdζ*dζdy_klm
                
                dudz = dudξ*dξdz_klm + dudη*dηdz_klm + dudζ*dζdz_klm
                dvdz = dvdξ*dξdz_klm + dvdη*dηdz_klm + dvdζ*dζdz_klm
                dwdz = dwdξ*dξdz_klm + dwdη*dηdz_klm + dwdζ*dζdz_klm
                
                # Transform potential temperature derivatives to physical coordinates
                dθdx = dθdξ*dξdx_klm + dθdη*dηdx_klm + dθdζ*dζdx_klm
                dθdy = dθdξ*dξdy_klm + dθdη*dηdy_klm + dθdζ*dζdy_klm
                dθdz = dθdξ*dξdz_klm + dθdη*dηdz_klm + dθdζ*dζdz_klm

                # Strain rate tensor (symmetric part of velocity gradient)
                S11 = dudx
                S22 = dvdy
                S33 = dwdz
                S12 = 0.5 * (dudy + dvdx)
                S13 = 0.5 * (dudz + dwdx)
                S23 = 0.5 * (dvdz + dwdy)

                # Rotation tensor (anti-symmetric part)
                Ω12 = 0.5 * (dudy - dvdx)
                Ω13 = 0.5 * (dudz - dwdx)
                Ω21 = -Ω12
                Ω23 = 0.5 * (dvdz - dwdy)
                Ω31 = -Ω13
                Ω32 = -Ω23
                
                # Strain rate magnitude
                Sij = sqrt(0.5 * (S11*S11 + S22*S22 + S33*S33) + (S12*S12 + S13*S13 + S23*S23))
                S2  = Sij*Sij
                
                # Filter width calculation
                Δ2      = Δ * Δ
                # @info Δ
                
                # Richardson number calculation for stability correction
                # Get reference potential temperature (local value)
                θ_ref = uprimitive[k,l,m,5]
                
                # Buoyancy frequency squared: N² = (g/θ) * dθ/dz
                # Note: assuming z is vertical (modify if different coordinate system)
                N2 = abs(θ_ref) > 1e-12 ? (g / θ_ref) * dθdz : 0.0
                
               
                # Richardson number: Ri = N²/S²
                Ri = (S2 > 1e-12 && N2 >= 0.0) ? N2 / S2 : 0.0
                
                # Stability function for Richardson correction
                # Various formulations exist; using a smooth transition
                f_Ri = if Ri >= Ri_crit
                    # Stable stratification suppresses turbulence
                    0.0
                elseif Ri >= 0.0
                    # Stable but sub-critical: reduce mixing
                    (1.0 - Ri/Ri_crit)^2
                else
                    # Unstable stratification: enhance mixing
                    sqrt(1.0 - 16.0*Ri)  # Ri is negative, so this increases mixing
                end
                
                # Apply Richardson stability correction to eddy viscosity
                # Base Smagorinsky eddy viscosity
                ν_t_base = (C_s * Δ)^2 * Sij
                
                # Richardson-corrected eddy viscosity
                ν_t = ν_t_base #* f_Ri
                
                # Compute scalar gradient for diffusion
                dqdξ = 0.0; dqdη = 0.0; dqdζ = 0.0
                
                @turbo for ii = 1:ngl
                    dqdξ += dψ[ii,k]*uprimitive[ii,l,m,ieq]
                    dqdη += dψ[ii,l]*uprimitive[k,ii,m,ieq]
                    dqdζ += dψ[ii,m]*uprimitive[k,l,ii,ieq]
                end
                
                # Transform scalar gradient to physical coordinates
                dqdx_phys = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
                dqdy_phys = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
                dqdz_phys = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
                
                # Determine effective diffusivity based on scalar type
                # TODO: Replace this logic with proper equation identification
                # Common orderings:
                # - Conservative: [ρ, ρu, ρv, ρw, ρE] or [ρ, ρu, ρv, ρw, ρE, ρθ]
                # - Primitive: [ρ, u, v, w, T] or [ρ, u, v, w, p, θ]
                if ieq == 5  # Assuming potential temperature equation is at index 5
                    # For temperature: use thermal diffusivity (ν_t / Pr_t)
                    effective_diffusivity = visc_coeffieq[ieq] * ν_t / Pr_t
                else
                    # For momentum equations: use momentum diffusivity
                    effective_diffusivity = visc_coeffieq[ieq] * ν_t
                end
                
                # Apply effective diffusivity to scalar gradients
                dqdx = effective_diffusivity * dqdx_phys
                dqdy = effective_diffusivity * dqdy_phys
                dqdz = effective_diffusivity * dqdz_phys
                
                # Transform back to computational coordinates for weak form
                ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz)*ωJac
                ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz)*ωJac
                ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz)*ωJac 
                
                
                # Distribute to element RHS arrays
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

function  _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeff, ω, mesh, basis, metrics, inputs, rhs_el, iel, ieq, QT::Exact, VT, SD::NSD_2D, ::FD)
    nothing
end

function  _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, 
                           uprimitive, visc_coeffieq, ω,
                           ngl, dψ, Je,
                           dξdx, dξdy, 
                           dηdx, dηdy, 
                           inputs,
                           rhs_el,
                           iel, ieq,
                           τ_f, wθ, lwall_model,
                           connijk,
                           coords,
                           poin_in_bdy_face, elem_to_face, bdy_face_type,
                           QT::Inexact, VT::AV, SD::NSD_2D, ::ContGal; Δ=1.0)
    
    nothing
    
end

function compute_vertical_derivative_q!(dqdz, q, iel, ngl, Je, dξdz, dηdz, dζdz, ω, dψ, ::NSD_3D)

    for k=1:ngl
        for j=1:ngl
            for i=1:ngl
                ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                
                dHdξ = 0.0
                dHdη = 0.0
                dHdζ = 0.0
                @turbo for m = 1:ngl
                    dHdξ += dψ[m,i]*q[m,j,k]
                    dHdη += dψ[m,j]*q[i,m,k]
                    dHdζ += dψ[m,k]*q[i,j,m]
                end
                dξdz_ij = dξdz[iel,i,j,k]
                dηdz_ij = dηdz[iel,i,j,k]
                dζdz_ij = dζdz[iel,i,j,k]
                
                dHdz = dHdξ*dξdz_ij + dHdη*dηdz_ij + dHdζ*dζdz_ij

                auxi = ωJac*dHdz
                dqdz[iel,i,j,k] += auxi
            end
        end
    end
end

function compute_vertical_derivative_q!(dqdz, q, iel, ngl, Je, dξdy, dηdy, ω, dψ, ::NSD_2D)
    for j=1:ngl
        for i=1:ngl
            ωJac = ω[i]*ω[j]*Je[iel,i,j]
            
            dHdξ = 0.0    
            dHdη = 0.0
            @turbo for m = 1:ngl
                dHdξ += dψ[m,i]*q[m,j]
                dHdη += dψ[m,j]*q[i,m]
            end
            dξdy_ij = dξdy[iel,i,j]      
            dηdy_ij = dηdy[iel,i,j]      
            
            dHdz = dHdξ*dξdy_ij + dHdη*dηdy_ij
            
            auxi = ωJac*dHdz
            dqdz[iel,i,j] += auxi
        end 
    end     
end  

function saturation_adjustment(uaux, qe, z, connijk, nelem, ngl, neqs, thermo_params)
    for iel=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    ip = connijk[iel,k,j,i]
                    @inbounds uaux[ip, 1:neqs] .= user_saturation_adjustment(@view(uaux[ip,:]), @view(qe[ip,:]), z[ip], thermo_params)
                end
            end
        end
    end
end
