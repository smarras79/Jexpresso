using Distributions
using StaticArrays

const PHYS_CONST = PhysicalConst{Float64}()

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

function micro2rhs!(rhs,dhldt,dqtdt,dqpdt,::NSD_2D)

    @view(rhs[:,:,:,4]) .= @view(rhs[:,:,:,4]) .- @view(dhldt[:,:,:])
    @view(rhs[:,:,:,5]) .= @view(rhs[:,:,:,5]) .+ @view(dqtdt[:,:,:])
    @view(rhs[:,:,:,6]) .= @view(rhs[:,:,:,6]) .+ @view(dqpdt[:,:,:])

end


function micro2rhs!(rhs,dhldt,dqtdt,dqpdt,::NSD_3D)

    @view(rhs[:,:,:,:,5]) .= @view(rhs[:,:,:,:,5]) .- @view(dhldt[:,:,:,:])
    @view(rhs[:,:,:,:,6]) .= @view(rhs[:,:,:,:,6]) .+ @view(dqtdt[:,:,:,:])
    @view(rhs[:,:,:,:,7]) .= @view(rhs[:,:,:,:,7]) .+ @view(dqpdt[:,:,:,:])

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
    # for @timers, do not delete
    timers = params.timers
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
            DSS_global_RHS!(@view(params.RHS[:,:]), params.g_dss_cache, params.neqs)

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
    Δt      = params.Δt
    # for @timers, do not delete
    timers  = params.timers

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
                    connijk_lag = params.mesh.connijk_lag, Je_lag = params.metrics_lag.Je, ladapt = inputs[:ladapt])
        else
            filter!(u, params, time, params.uaux, params.mesh.connijk, params.metrics.Je, SD, params.SOL_VARS_TYPE; ladapt = inputs[:ladapt])
        end
    end

    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    if inputs[:ladapt] == true
        conformity4ncf_q!(params.uaux, params.rhs_el_tmp, @view(params.utmp[:,1:neqs]), params.vaux, 
                            params.g_dss_cache,
                            params.mesh.SD, 
                            params.QT, params.mesh.connijk,
                            params.mesh, params.Minv, 
                            params.metrics.Je, params.ω, params.AD, 
                            params.neqs,
                            params.q_el, params.q_el_pro,
                            params.cache_ghost_p, params.q_ghost_p,
                            params.cache_ghost_c, params.q_ghost_c,
                            params.interp)
    end
    
    resetbdyfluxToZero!(params)
    apply_boundary_conditions_dirichlet!(u, params.uaux, time, params.qp.qe,
                                         params.mesh.coords,
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
            #=@time if (SD == NSD_3D())
                params.rhs_el[:,:,:,:,5] .-= params.mp.dhldt
                params.rhs_el[:,:,:,:,6] .+= params.mp.dqtdt
                params.rhs_el[:,:,:,:,7] .+= params.mp.dqpdt
            else
                @time @view(params.rhs_el[:,:,:,4]) .-= params.mp.dhldt[:,:,:]
                @time @view(params.rhs_el[:,:,:,5]) .+= @view(params.mp.dqtdt[:,:,:])
                @time @view(params.rhs_el[:,:,:,6]) .= @view(params.rhs_el[:,:,:,6]) .+ @view(params.mp.dqpdt[:,:,:])
            end=#
            micro2rhs!(params.rhs_el,params.mp.dhldt, params.mp.dqtdt, params.mp.dqpdt, SD)
        end
        uaux2u!(u, params.uaux, params.neqs, params.mesh.npoin)
    end
    
    if(params.inputs[:lsaturation])
        saturation_adjustment(params.uaux, params.qp.qe, params.mesh.z, params.mesh.connijk, params.mesh.nelem, params.mesh.ngl, neqs, params.thermo_params)
        uaux2u!(u, params.uaux, params.neqs, params.mesh.npoin)
    end
    
    inviscid_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, params.mesh.coords, lsource, 
                     params.mp.S_micro, params.mp.qn, params.mp.flux_lw, params.mp.flux_sw, SD)
    
    if inputs[:ladapt] == true
        DSS_nc_gather_rhs!(params.RHS, SD, QT, params.rhs_el,
                           params.mesh.non_conforming_facets,
                           params.mesh.non_conforming_facets_parents_ghost, params.cache_ghost_p,
                           params.q_el, params.q_el_pro, params.q_ghost_p,
                           params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPc_list_pg,
                           params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost, params.mesh.pgip_local, ngl-1, neqs, params.interp)

    end
    DSS_rhs!(params.RHS, params.rhs_el, params.mesh.connijk, nelem, ngl, neqs, SD, AD)

    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)
        
        resetRHSToZero_viscous!(params, SD)
        
        viscous_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, SD)
        
        if inputs[:ladapt] == true
            DSS_nc_gather_rhs!(params.RHS_visc, SD, QT, params.rhs_diff_el,
                           params.mesh.non_conforming_facets,
                           params.mesh.non_conforming_facets_parents_ghost, params.cache_ghost_p,
                           params.q_el, params.q_el_pro, params.q_ghost_p,
                           params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPc_list_pg,
                           params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost, params.mesh.pgip_local, ngl-1, neqs, params.interp)

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
    
    DSS_global_RHS!(@view(params.RHS[:,:]), params.g_dss_cache, params.neqs)
    
    #if (rem(time, Δt) == 0 && time > 0.0)
    if (time > 0.0)
        params.qp.qnm1 .= params.qp.qnm2
        params.qp.qnm2 .= params.uaux
    end
    
    for ieq=1:neqs
        divide_by_mass_matrix!(@view(params.RHS[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
        
        if inputs[:ladapt] == true
            
            DSS_nc_scatter_rhs!(@view(params.RHS[:,ieq]), SD, QT,
                                params.mesh.non_conforming_facets,
                                params.mesh.non_conforming_facets_children_ghost, params.cache_ghost_c,
                                params.q_el, params.q_el_pro, params.q_ghost_c,
                                params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPp_list_cg,
                                params.mesh.gip2ip, params.mesh.cgip_local, ngl-1, params.interp)
        end
    end
end

function inviscid_rhs_el!(u, params, connijk, qe, coords, lsource, S_micro_vec, qn_vec, flux_lw_vec, flux_sw_vec, SD::NSD_1D)
    
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)

    ngl   = params.mesh.ngl
    npoin = params.mesh.npoin
    nelem = params.mesh.nelem
    neqs  = params.neqs
    
    xmin = params.xmin; xmax = params.xmax; ymax = params.ymax
    
    for iel=1:nelem   
        for i=1:ngl
            ip = connijk[iel,i,1]
            
            user_primitives!(@view(params.uaux[ip,:]),
                             @view(qe[ip,:]),
                             @view(params.uprimitive[i,:]),
                             params.SOL_VARS_TYPE)

            user_flux!(@view(params.F[i,:]),
                       @view(params.G[i,:]), SD,
                       @view(params.uaux[ip,:]),
                       @view(qe[ip,:]),
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip)
            
            if lsource
                user_source!(@view(params.S[i,:]),
                             @view(params.uaux[ip,:]),
                             @view(qe[ip,:]),
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE;
                             neqs=params.neqs, x=coords[ip,1], y=0.0, xmax=xmax,xmin=xmin)
            end
        end

        _expansion_inviscid!(u, params.neqs, ngl,
                             params.basis.dψ, params.ω,
                             params.F, params.S,
                             params.rhs_el,
                             iel, params.CL, params.QT, SD, params.AD)
        
    end
end

function inviscid_rhs_el!(u, params,
                          connijk::Array{Int64,4},
                          qe::Matrix{Float64},
                          coords, 
                          lsource, S_micro_vec, qn_vec, flux_lw_vec,
                          flux_sw_vec, SD::NSD_2D)
    
    PhysConst = PhysicalConst{Float64}()

    u_element_wise = zeros(params.mesh.ngl, params.mesh.ngl, params.neqs)
    
    xmin = params.xmin; xmax = params.xmax; ymax = params.ymax
    for iel = 1:params.mesh.nelem
        
        
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
                    S_micro::Float64 = @inbounds S_micro_vec[ip]
                    flux_lw::Float64 = @inbounds flux_lw_vec[ip]
                    flux_sw::Float64 = @inbounds flux_sw_vec[ip]
                    qn::Float64 = @inbounds qn_vec[ip]
                    add_micro_precip_sources!(@view(params.S[i,j,:]),
                                              @view(params.uaux[ip,:]),
                                              @view(qe[ip,:]),
                                              S_micro, qn, flux_lw, flux_sw, PHYS_CONST,
                                              SD, params.SOL_VARS_TYPE)
                end
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

function inviscid_rhs_el!(u, params, connijk, qe, coords, lsource, S_micro_vec, qn_vec, flux_lw_vec, flux_sw_vec, SD::NSD_3D)
   
    PhysConst = PhysicalConst{Float64}()
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
                    S_micro::Float64 = @inbounds S_micro_vec[ip]
                    flux_lw::Float64 = @inbounds flux_lw_vec[ip]
                    flux_sw::Float64 = @inbounds flux_sw_vec[ip]
                    qn::Float64 = @inbounds qn_vec[ip]
                    add_micro_precip_sources!(@view(params.S[i,j,k,:]),
                                                @view(params.uaux[ip,:]),
                                                @view(qe[ip,:]),
                                                S_micro, qn, flux_lw, flux_sw, PhysConst,
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
                             params.mesh.connijk,
                             params.mesh.coords,
                             params.mesh.poin_in_bdy_face,
                             params.mesh.elem_to_face,
                             params.mesh.bdy_face_type,
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
                             iel, ieq,
                             params.QT, params.VT, SD, params.AD; Δ=Δ)
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
                             params.inputs, params.rhs_el, iel, ieq, params.mesh.connijk,
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
                              connijk,
                              coords,
                              poin_in_bdy_face,
                              elem_to_face,
                              bdy_face_type,
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
# viscous RHS 2D
#
function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, ieq,
                          QT::Inexact, VT, SD::NSD_2D, ::ContGal; Δ=1.0, vargs...)

    PhysConst = PhysicalConst{Float32}()
    Sc_t      = PhysConst.Sc_t
    Δ2        = Δ^2

    # Determine if this is a momentum equation
    is_u_momentum  = (ieq == 2)
    is_v_momentum  = (ieq == 3)
    is_temperature = (ieq == 4)
    
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

            #∇⋅u
            div_u = dudx + dvdy

            if is_u_momentum
                # USE EFFECTIVE VISCOSITY
                effective_viscosity =  SGS_diffusion(visc_coeffieq, ieq,
                                                     uprimitiveieq[k,l,1],
                                                     dudx, dvdy, dudy, dvdx,
                                                     PhysConst, Δ2,
                                                     inputs, 
                                                     VT, SD)
                
                τ_xx = 2.0 * effective_viscosity * dudx - (2.0/3.0) * effective_viscosity * div_u
                τ_xy = effective_viscosity * (dudy + dvdx)
                flux_x = τ_xx
                flux_y = τ_xy

                
            elseif is_v_momentum
                # USE EFFECTIVE VISCOSITY
                effective_viscosity =  SGS_diffusion(visc_coeffieq, ieq,
                                                     uprimitiveieq[k,l,1],
                                                     dudx, dvdy, dudy, dvdx,
                                                     PhysConst, Δ2,
                                                     inputs, 
                                                     VT, SD)
                
                τ_xy = effective_viscosity * (dudy + dvdx)
                τ_yy = 2.0 * effective_viscosity * dvdy - (2.0/3.0) * effective_viscosity * div_u
                flux_x = τ_xy
                flux_y = τ_yy

                
            elseif is_temperature
                # USE EFFECTIVE DIFFUSIVITY
                effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                      uprimitiveieq[k,l,1],
                                                      dudx, dvdy, dudy, dvdx,
                                                      PhysConst, Δ2,
                                                      inputs, 
                                                      VT, SD)
                
                # Compute temperature gradient
                dθdξ = 0.0; dθdη = 0.0
                @turbo for ii = 1:ngl
                    dθdξ += dψ[ii,k]*uprimitiveieq[ii,l,ieq]
                    dθdη += dψ[ii,l]*uprimitiveieq[k,ii,ieq]
                end
                
                dθdx = dθdξ*dξdx_kl + dθdη*dηdx_kl
                dθdy = dθdξ*dξdy_kl + dθdη*dηdy_kl
                
                flux_x = effective_diffusivity * dθdx
                flux_y = effective_diffusivity * dθdy
                
            else
                # Other scalars (use appropriate Schmidt number)
                # USE EFFECTIVE DIFFUSIVITY
                effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                      uprimitiveieq[k,l,1],
                                                      dudx, dvdy, dudy, dvdx,
                                                      PhysConst, Δ2,
                                                      inputs, 
                                                      VT, SD)
                
                # Compute temperature gradient
                dqdξ = 0.0; dqdη = 0.0
                @turbo for ii = 1:ngl
                    dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,ieq]
                    dqdη += dψ[ii,l]*uprimitiveieq[k,ii,ieq]
                end
                
                dqdx = dqdξ*dξdx_kl + dqdη*dηdx_kl
                dqdy = dqdξ*dξdy_kl + dqdη*dηdy_kl
                
                flux_x = effective_diffusivity * dqdx
                flux_y = effective_diffusivity * dqdy
            end

            # ===== Weak form assembly (same for all) =====
            ∇ξ_flux_kl = (dξdx_kl*flux_x + dξdy_kl*flux_y)*ωJac
            ∇η_flux_kl = (dηdx_kl*flux_x + dηdy_kl*flux_y)*ωJac
            
            @turbo for i = 1:ngl
                dhdξ_ik = dψ[i,k]
                dhdη_il = dψ[i,l]
                
                rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * ∇ξ_flux_kl
                rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * ∇η_flux_kl
            end
        end  
    end
end


function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy, dξdz,
                          dηdx, dηdy, dηdz,
                          dζdx, dζdy, dζdz,
                          inputs,
                          rhs_el,
                          iel, ieq,
                          connijk,
                          coords, 
                          poin_in_bdy_face, elem_to_face, bdy_face_type,
                          QT::Inexact, VT::AV, SD::NSD_3D, ::ContGal; Δ=1.0)

    PhysConst = PhysicalConst{Float32}()
    
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



function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy, dξdz,
                          dηdx, dηdy, dηdz,
                          dζdx, dζdy, dζdz,
                          inputs, rhs_el,
                          iel, ieq, connijk,
                          coords, 
                          poin_in_bdy_face, elem_to_face, bdy_face_type,
                          QT::Inexact, VT, SD::NSD_3D, ::ContGal; Δ=1.0)

    
    PhysConst = PhysicalConst{Float32}()
    Δ2        = Δ^2

    # Determine equation type (indices shifted for 3D)
    is_u_momentum  = (ieq == 2)
    is_v_momentum  = (ieq == 3)
    is_w_momentum  = (ieq == 4)
    is_temperature = (ieq == 5)
    
    for m = 1:ngl      # ADDED: third loop
        for l = 1:ngl
            for k = 1:ngl
                ωJac = ω[k]*ω[l]*ω[m]*Je[iel,k,l,m]

                # ===== Compute all velocity gradients =====
                # u-velocity gradients
                dudξ = 0.0; dudη = 0.0; dudζ = 0.0
                @turbo for ii = 1:ngl
                    dudξ += dψ[ii,k]*uprimitiveieq[ii,l,m,2]
                    dudη += dψ[ii,l]*uprimitiveieq[k,ii,m,2]
                    dudζ += dψ[ii,m]*uprimitiveieq[k,l,ii,2]
                end
                
                # v-velocity gradients
                dvdξ = 0.0; dvdη = 0.0; dvdζ = 0.0
                @turbo for ii = 1:ngl
                    dvdξ += dψ[ii,k]*uprimitiveieq[ii,l,m,3]
                    dvdη += dψ[ii,l]*uprimitiveieq[k,ii,m,3]
                    dvdζ += dψ[ii,m]*uprimitiveieq[k,l,ii,3]
                end
                
                # w-velocity gradients (NEW)
                dwdξ = 0.0; dwdη = 0.0; dwdζ = 0.0
                @turbo for ii = 1:ngl
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

                # Transform to physical coordinates
                # u-velocity
                dudx = dudξ*dξdx_klm + dudη*dηdx_klm + dudζ*dζdx_klm 
                dudy = dudξ*dξdy_klm + dudη*dηdy_klm + dudζ*dζdy_klm
                dudz = dudξ*dξdz_klm + dudη*dηdz_klm + dudζ*dζdz_klm
                
                # v-velocity
                dvdx = dvdξ*dξdx_klm + dvdη*dηdx_klm + dvdζ*dζdx_klm
                dvdy = dvdξ*dξdy_klm + dvdη*dηdy_klm + dvdζ*dζdy_klm
                dvdz = dvdξ*dξdz_klm + dvdη*dηdz_klm + dvdζ*dζdz_klm
                
                # w-velocity (NEW)
                dwdx = dwdξ*dξdx_klm + dwdη*dηdx_klm + dwdζ*dζdx_klm
                dwdy = dwdξ*dξdy_klm + dwdη*dηdy_klm + dwdζ*dζdy_klm
                dwdz = dwdξ*dξdz_klm + dwdη*dηdz_klm + dwdζ*dζdz_klm

                # Velocity divergence
                div_u = dudx + dvdy + dwdz

                if is_u_momentum
                    # USE EFFECTIVE VISCOSITY
                    effective_viscosity = SGS_diffusion(visc_coeffieq, ieq,
                                                        uprimitiveieq[k,l,m,1],
                                                        dudx, dvdy, dwdz,      
                                                        dudy, dvdx,            
                                                        dudz, dwdx,            
                                                        dvdz, dwdy,
                                                        0.0,
                                                        0.0,
                                                        PhysConst, Δ2,
                                                        inputs, 
                                                        VT, SD)
                    
                    # Stress tensor for u-momentum
                    τ_xx = 2.0 * effective_viscosity * dudx - (2.0/3.0) * effective_viscosity * div_u
                    τ_xy = effective_viscosity * (dudy + dvdx)
                    τ_xz = effective_viscosity * (dudz + dwdx)
                    
                    flux_x = τ_xx
                    flux_y = τ_xy
                    flux_z = τ_xz

                    
                elseif is_v_momentum
                    # USE EFFECTIVE VISCOSITY
                    effective_viscosity = SGS_diffusion(visc_coeffieq, ieq,
                                                        uprimitiveieq[k,l,m,1],
                                                        dudx, dvdy, dwdz,      
                                                        dudy, dvdx,            
                                                        dudz, dwdx,            
                                                        dvdz, dwdy, 
                                                        0.0,
                                                        0.0,           
                                                        PhysConst, Δ2,
                                                        inputs, 
                                                        VT, SD)
                    
                    # Stress tensor for v-momentum
                    τ_xy = effective_viscosity * (dudy + dvdx)
                    τ_yy = 2.0 * effective_viscosity * dvdy - (2.0/3.0) * effective_viscosity * div_u
                    τ_yz = effective_viscosity * (dvdz + dwdy)
                    
                    flux_x = τ_xy
                    flux_y = τ_yy
                    flux_z = τ_yz

                    
                elseif is_w_momentum  # NEW BLOCK
                    # USE EFFECTIVE VISCOSITY
                    effective_viscosity = SGS_diffusion(visc_coeffieq, ieq,
                                                        uprimitiveieq[k,l,m,1],
                                                        dudx, dvdy, dwdz,
                                                        dudy, dvdx,
                                                        dudz, dwdx,
                                                        dvdz, dwdy,
                                                        0.0,
                                                        0.0,
                                                        PhysConst, Δ2,
                                                        inputs, 
                                                        VT, SD)
                    
                    # Stress tensor for w-momentum
                    τ_xz = effective_viscosity * (dudz + dwdx)
                    τ_yz = effective_viscosity * (dvdz + dwdy)
                    τ_zz = 2.0 * effective_viscosity * dwdz - (2.0/3.0) * effective_viscosity * div_u
                    
                    flux_x = τ_xz
                    flux_y = τ_yz
                    flux_z = τ_zz

                    
                elseif is_temperature
                   
                    # Compute temperature gradient
                    dθdξ = 0.0; dθdη = 0.0; dθdζ = 0.0
                    @turbo for ii = 1:ngl
                        dθdξ += dψ[ii,k]*uprimitiveieq[ii,l,m,ieq]
                        dθdη += dψ[ii,l]*uprimitiveieq[k,ii,m,ieq]
                        dθdζ += dψ[ii,m]*uprimitiveieq[k,l,ii,ieq]
                    end
                    
                    # Transform to physical coordinates
                    dθdx = dθdξ*dξdx_klm + dθdη*dηdx_klm + dθdζ*dζdx_klm
                    dθdy = dθdξ*dξdy_klm + dθdη*dηdy_klm + dθdζ*dζdy_klm
                    dθdz = dθdξ*dξdz_klm + dθdη*dηdz_klm + dθdζ*dζdz_klm

                    if inputs[:energy_equation] == "theta" && inputs[:lrichardson]
                        θ_ref = uprimitiveieq[k,l,m,5]  # Local temperature
                    else
                        θ_ref = 1.0  # Dummy value (not used when lrichardson=false)
                    end
                    
                    # USE EFFECTIVE DIFFUSIVITY
                    effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                          uprimitiveieq[k,l,m,1],
                                                          dudx, dvdy, dwdz,      
                                                          dudy, dvdx,            
                                                          dudz, dwdx,            
                                                          dvdz, dwdy,
                                                          θ_ref,
                                                          dθdz,
                                                          PhysConst, Δ2,
                                                          inputs, 
                                                          VT, SD)
                    
                    
                    flux_x = effective_diffusivity * dθdx
                    flux_y = effective_diffusivity * dθdy
                    flux_z = effective_diffusivity * dθdz
                    
                else
                    # Other scalars (use appropriate Schmidt number)
                    # USE EFFECTIVE DIFFUSIVITY
                    effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                          uprimitiveieq[k,l,m,1],
                                                          dudx, dvdy, dwdz,      
                                                          dudy, dvdx,            
                                                          dudz, dwdx,            
                                                          dvdz, dwdy,
                                                          0.0,
                                                          0.0,
                                                          PhysConst, Δ2,
                                                          inputs, 
                                                          VT, SD)
                    
                    # Compute scalar gradient
                    dqdξ = 0.0; dqdη = 0.0; dqdζ = 0.0
                    @turbo for ii = 1:ngl
                        dqdξ += dψ[ii,k]*uprimitiveieq[ii,l,m,ieq]
                        dqdη += dψ[ii,l]*uprimitiveieq[k,ii,m,ieq]
                        dqdζ += dψ[ii,m]*uprimitiveieq[k,l,ii,ieq]
                    end
                    
                    # Transform to physical coordinates
                    dqdx = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
                    dqdy = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
                    dqdz = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
                    
                    flux_x = effective_diffusivity * dqdx
                    flux_y = effective_diffusivity * dqdy
                    flux_z = effective_diffusivity * dqdz
                end

                # ===== Weak form assembly (3D) =====
                ∇ξ_flux_klm = (dξdx_klm*flux_x + dξdy_klm*flux_y + dξdz_klm*flux_z)*ωJac
                ∇η_flux_klm = (dηdx_klm*flux_x + dηdy_klm*flux_y + dηdz_klm*flux_z)*ωJac
                ∇ζ_flux_klm = (dζdx_klm*flux_x + dζdy_klm*flux_y + dζdz_klm*flux_z)*ωJac
                
                @turbo for i = 1:ngl
                    dhdξ_ik = dψ[i,k]
                    dhdη_il = dψ[i,l]
                    dhdζ_im = dψ[i,m]
                    
                    rhs_diffξ_el[iel,i,l,m,ieq] -= dhdξ_ik * ∇ξ_flux_klm
                    rhs_diffη_el[iel,k,i,m,ieq] -= dhdη_il * ∇η_flux_klm
                    rhs_diffζ_el[iel,k,l,i,ieq] -= dhdζ_im * ∇ζ_flux_klm
                end
            end
        end  
    end
end
##


function  _expansion_visc_old!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                           uprimitiveieq, visc_coeffieq, ω,
                           ngl, dψ, Je,
                           dξdx, dξdy, dξdz,
                           dηdx, dηdy, dηdz,
                           dζdx, dζdy, dζdz,
                           inputs,
                           rhs_el,
                           iel, ieq,
                           connijk,
                           coords, 
                           poin_in_bdy_face, elem_to_face, bdy_face_type,
                           QT::Inexact, VT, SD::NSD_3D, ::ContGal; Δ=1.0, vargs...)
    
    PhysConst  = PhysicalConst{Float32}()
    Δ2         = Δ * Δ
    for m = 1:ngl
        for l = 1:ngl
            for k = 1:ngl
                Je_klm = Je[iel,k,l,m]
                ωJac = ω[k]*ω[l]*ω[m]*Je_klm
                
                # Velocity gradients in computational space
                dudξ = 0.0; dudη = 0.0; dudζ = 0.0
                dvdξ = 0.0; dvdη = 0.0; dvdζ = 0.0
                dwdξ = 0.0; dwdη = 0.0; dwdζ = 0.0
                dθdξ = 0.0; dθdη = 0.0; dθdζ = 0.0
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
                    
                    # Potential temperature gradients
                    dθdξ += dψ[ii,k]*uprimitiveieq[ii,l,m,5]
                    dθdη += dψ[ii,l]*uprimitiveieq[k,ii,m,5]
                    dθdζ += dψ[ii,m]*uprimitiveieq[k,l,ii,5]
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

                # Transform potential temperature derivatives to physical coordinates
                dθdz  = dθdξ*dξdz_klm + dθdη*dηdz_klm + dθdζ*dζdz_klm

                θ_ref = uprimitiveieq[k,l,m,5]
                
                effective_diffusivity = SGS_diffusivity(visc_coeffieq, ieq,
                                                        uprimitiveieq[k,l,m,1], 
                                                        u11, u12, u13,
                                                        u21, u22, u23,
                                                        u31, u32, u33,
                                                        θ_ref,
                                                        dθdz, 
                                                        PhysConst, Δ2,
                                                        inputs, 
                                                        VT, SD)
                
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

function compute_vertical_derivative_q!(dqdz::Array{Float64,4}, q::Array{Float64,4}, iel::Int64, ngl::Int64, Je::Array{Float64,4}, 
        dξdz::Array{Float64,4}, dηdz::Array{Float64,4}, dζdz::Array{Float64,4}, ω::Vector{Float64}, dψ::Matrix{Float64}, ::NSD_3D)
    
    local ωJac::Float64
    local dHdξ::Float64
    local dHdη::Float64
    local dHdζ::Float64
    local dξdz_ij::Float64
    local dηdz_ij::Float64
    local dζdz_ij::Float64
    local dHdz::Float64
    local auxi::Float64

    for k=1:ngl
        for j=1:ngl
            for i=1:ngl
                @inbounds ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                
                dHdξ = 0.0
                dHdη = 0.0
                dHdζ = 0.0
                @turbo for m = 1:ngl
                    dHdξ += dψ[m,i]*q[m,j,k,1]
                    dHdη += dψ[m,j]*q[i,m,k,1]
                    dHdζ += dψ[m,k]*q[i,j,m,1]
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

function compute_vertical_derivative_q!(dqdz, q, iel::Int64, ngl::Int64, Je, dξdy, dηdy, ω, dψ, ::NSD_2D)
    for j=1:ngl
        for i=1:ngl
            ωJac = ω[i]*ω[j]*Je[iel,i,j]
            
            dHdξ = 0.0    
            dHdη = 0.0
            @turbo for m = 1:ngl
                dHdξ += dψ[m,i]*q[m,j,1]
                dHdη += dψ[m,j]*q[i,m,1]
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

@inline function flux_ranocha(u_ll, u_rr, orientation::Integer,
                              equations)
    # Unpack left and right state
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
    params.uprimitive[i,:]

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
    # in exact arithmetic since
    #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
    #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
    v1_avg = 0.5 * (v1_ll + v1_rr)
    p_avg = 0.5 * (p_ll + p_rr)
    velocity_square_avg = 0.5 * (v1_ll * v1_rr)
    
    # Calculate fluxes
    # Ignore orientation since it is always "1" in 1D
    f1 = rho_mean * v1_avg
    f2 = f1 * v1_avg + p_avg
    f3 = f1 * (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
         0.5 * (p_ll * v1_rr + p_rr * v1_ll)

    return SVector(f1, f2, f3)
end
