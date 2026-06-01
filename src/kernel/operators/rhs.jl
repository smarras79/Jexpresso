using Distributions
using TrixiBase
using StaticArrays
using UnPack

# TODO - This hardcoded 4 is terrible!
# Even worse. The indices are assuming that the order is the opposite of what high performance requires
@inline @inbounds function get_node_vars_4(u, index1)
    # SVector(ntuple(@inline(v->u[indices..., v]), 4))
    SVector(u[index1, 1], u[index1, 2], u[index1, 3], u[index1, 4])
end

@inline @inbounds function get_node_vars_4(u, index1, index2)
    # SVector(ntuple(@inline(v->u[indices..., v]), 4))
    SVector(u[index1, index2, 1], u[index1, index2, 2], u[index1, index2, 3], u[index1, index2, 4])
end

@inline @inbounds function set_node_vars_4!(u, values, index1)
    u[index1, 1] = values[1]
    u[index1, 2] = values[2]
    u[index1, 3] = values[3]
    u[index1, 4] = values[4]
end

@inline @inbounds function set_node_vars_4!(u, values, index1, index2)
    u[index1, index2, 1] = values[1]
    u[index1, index2, 2] = values[2]
    u[index1, index2, 3] = values[3]
    u[index1, index2, 4] = values[4]
end

const PHYS_CONST = PhysicalConst{Float64}()
const MicroConst = MicrophysicalConst{Float64}()

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
@trixi_timeit timer() "rhs" begin
    # backend = params.inputs[:backend]
    backend = CPU()
    # for @timers, do not delete
    timers = params.timers
    if (backend == CPU())
        _build_rhs!(params.RHS, u, params, time)

        if (params.laguerre)
            build_rhs_laguerre!(params.RHS_lag, u, params, time)
            params.RHS .+= params.RHS_lag
        end

        RHStoDU!(du, params.RHS, params.neqs, params.mesh.npoin)
    else
        if (params.SOL_VARS_TYPE == PERT())
            lpert = true
        else
            lpert = false
        end

        if (params.SD == NSD_1D())
            params.RHS .= TFloat(0.0)
            k1 = utouaux_gpu!(backend)
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
               workgroupsize = (params.neqs))

            k = _build_rhs_gpu_v0!(backend,(Int64(params.mesh.ngl)))
            k(params.RHS, u, params.uaux, params.qp.qe, params.mesh.x, TFloat(time),
              params.mesh.connijk , params.basis.dψ, params.ω, params.Minv,
              params.flux_gpu, params.source_gpu,
              PHYS_CONST, params.xmax, params.xmin, params.mesh.ngl, params.neqs,
              lpert, inputs[:lperiodic_1d], params.mesh.npoin_linear, params.mesh.npoin;
              ndrange = params.mesh.nelem*params.mesh.ngl,workgroupsize = params.mesh.ngl)

            if (params.laguerre)
                params.RHS_lag .= TFloat(0.0)
                k = _build_rhs_gpu_v0!(backend,(Int64(params.mesh.ngr)))
                k(params.RHS, u, params.uaux, params.qp.qe, params.mesh.x, TFloat(time),
                  params.mesh.connijk_lag , params.basis_lag.dψ, params.ω_lag, params.Minv,
                  params.flux_lag_gpu, params.source_lag_gpu,
                  PHYS_CONST, params.xmax, params.xmin, params.mesh.ngr, params.neqs,
                  lpert, inputs[:lperiodic_1d], params.mesh.npoin_linear, params.mesh.npoin;
                  ndrange = params.mesh.nelem_semi_inf*params.mesh.ngr,workgroupsize = params.mesh.ngr)

                @inbounds  params.RHS .+= params.RHS_lag
            end
            k1 = RHStodu_gpu!(backend)
            k1(params.RHS,du,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
               workgroupsize = (params.mesh.ngl,params.neqs))
        elseif (params.SD == NSD_3D())

            params.RHS .= TFloat(0.0)

            k1 = utouaux_gpu!(backend)
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
               workgroupsize = (params.neqs))

            if (params.inputs[:lfilter])
                params.B .= TFloat(0.0)
                kf = filter_gpu_3d!(backend,(Int64(params.mesh.ngl), Int64(params.mesh.ngl), Int64(params.mesh.ngl)))
                kf(@view(params.uaux[:,:]), params.qp.qe, params.B, params.fx, params.fy_t, params.fz_t,
                   params.metrics.Je, params.ω, params.ω, params.ω, params.mesh.connijk, params.Minv,
                   params.mesh.ngl, params.mesh.ngl, params.mesh.ngl, params.neqs, lpert;
                   ndrange = (params.mesh.nelem * params.mesh.ngl, params.mesh.ngl, params.mesh.ngl),
                   workgroupsize = (params.mesh.ngl, params.mesh.ngl, params.mesh.ngl))

                KernelAbstractions.synchronize(backend)
                if (lpert)
                    params.uaux[:,1:params.neqs] .= params.B
                else
                    params.uaux .= params.B .+ params.qp.qe
                end
                kf = uauxtou_gpu!(backend)
                kf(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
                   workgroupsize = (params.mesh.ngl,params.neqs))

                KernelAbstractions.synchronize(backend)
            end

            k = apply_boundary_conditions_gpu_3D!(backend)
            k(@view(params.uaux[:,:]), @view(u[:]), params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z,
              TFloat(time),params.metrics.nx,params.metrics.ny, params.metrics.nz,
              params.mesh.poin_in_bdy_face,params.qbdy_gpu,params.mesh.ngl,TInt(params.neqs),
              params.mesh.npoin, lpert;
              ndrange = (params.mesh.nfaces_bdy*params.mesh.ngl,params.mesh.ngl),
              workgroupsize = (params.mesh.ngl,params.mesh.ngl))

            KernelAbstractions.synchronize(backend)

            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs), workgroupsize = (params.neqs))

            if (inputs[:lmoist])
                k_moist = do_micro_physics_gpu_3D!(backend)
                k_moist(@view(params.uaux[:,:]), params.qp.qe, params.mp.Tabs, params.mp.qn, params.mp.qi, params.mp.qc,
                        params.mp.qr, params.mp.qs, params.mp.qg, params.mp.Pr, params.mp.Ps, params.mp.Pg,
                        params.mp.S_micro, PHYS_CONST, MicroConst, lpert, params.neqs, params.mesh.npoin, params.mesh.z,
                        params.adjusted, params.Pm; ndrange = (params.mesh.npoin))

                k_precip = _build_precipitation_rhs_gpu_3D_v0!(backend, (Int64(params.mesh.ngl),
                                                                         Int64(params.mesh.ngl),
                                                                         Int64(params.mesh.ngl)))

                k_precip(params.RHS, @view(params.uaux[:,:]), params.qp.qe,
                         params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk,
                         params.metrics.dξdz, params.metrics.dηdz, params.metrics.dζdz, params.metrics.Je,
                         params.basis.dψ, params.ω, params.Minv, params.flux_micro, params.source_micro,
                         params.mesh.ngl, TInt(params.neqs), PHYS_CONST, params.mesh.xmax, params.mesh.xmin,
                         params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert,
                         params.mp.Pr, params.mp.Ps, params.mp.Pg, params.mp.qi, params.mp.qn, params.mp.Tabs,
                         params.mp.S_micro, MicroConst;
                         ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl),
                         workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
            end
            KernelAbstractions.synchronize(backend)
            k = _build_rhs_gpu_3D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
            k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z,
              params.mesh.connijk, params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx,
              params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz,
              params.metrics.Je,
              params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu,
              params.mesh.ngl, TInt(params.neqs), PHYS_CONST,
              params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert;
              ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl),
              workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
            if (params.inputs[:case] != "bomex")
                k = _build_rhs_gpu_3D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z,
                  params.mesh.connijk, params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx,
                  params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz,
                  params.metrics.Je,
                  params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu,
                  params.mesh.ngl, TInt(params.neqs), PHYS_CONST,
                  params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert;
                  ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl),
                  workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
            else
                k = _build_rhs_gpu_3D_v1!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.z,
                  params.mesh.connijk, params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz, params.metrics.dηdx,
                  params.metrics.dηdy, params.metrics.dηdz, params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz,
                  params.metrics.Je,
                  params.basis.dψ, params.ω, params.Minv, params.flux_gpu, params.source_gpu,
                  params.mesh.ngl, TInt(params.neqs), PHYS_CONST,
                  params.thermo_params,
                  params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, params.mesh.zmax, params.mesh.zmin, lpert;
                  ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl),
                  workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))
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
                    k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el, params.rhs_diffζ_el,
                      params.uaux, params.qp.qe, params.source_gpu,
                      params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk,
                      params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz,
                      params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz,
                      params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz,
                      params.metrics.Je,
                      params.basis.dψ, params.ω, params.Minv,
                      params.visc_coeff,
                      params.mesh.ngl, TInt(params.neqs), PHYS_CONST, lpert;
                      ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl),
                      workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))

                elseif params.VT == SMAG()
                    k = _build_rhs_diff_gpu_3D_smag!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
                    k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el, params.rhs_diffζ_el,
                      params.uaux, params.qp.qe, params.source_gpu,
                      params.mesh.x, params.mesh.y, params.mesh.z, params.mesh.connijk,
                      params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz,
                      params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz,
                      params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz,
                      params.metrics.Je, params.basis.dψ, params.ω, params.Minv,
                      params.visc_coeff,
                      params.mesh.ngl, TInt(params.neqs), params.mesh.Δeffective_s, PHYS_CONST, lpert;
                      ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl),
                      workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))

                end
                KernelAbstractions.synchronize(backend)
                if (params.inputs[:case] == "bomex")
                    # param_set = TP.ThermodynamicsParameters(TFloat)
                    k_sa = saturation_adjustment_gpu_3D!(backend)
                    k_sa(params.uaux, params.qp.qe, params.mesh.z, params.mesh.connijk, TInt(params.neqs), params.thermo_params, lpert;
                         ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl,params.mesh.ngl),
                         workgroupsize = (params.mesh.ngl,params.mesh.ngl,params.mesh.ngl))

                    KernelAbstractions.synchronize(backend)

                    kf = uauxtou_gpu!(backend)
                    kf(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
                       workgroupsize = (params.mesh.ngl,params.neqs))
                    KernelAbstractions.synchronize(backend)
                end

                @inbounds params.RHS .+= params.RHS_visc
            end
            KernelAbstractions.synchronize(backend)

            k1 = RHStodu_gpu!(backend)
            k1(params.RHS,du,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
               workgroupsize = (params.mesh.ngl,params.neqs))

        elseif (params.SD == NSD_2D())
            params.RHS .= TFloat(0.0)

            k1 = utouaux_gpu!(backend)
            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
               workgroupsize = (params.mesh.ngl, params.neqs))

            if (params.inputs[:lfilter])
                params.B .= TFloat(0.0)
                kf = filter_gpu_2d!(backend,(Int64(params.mesh.ngl), Int64(params.mesh.ngl)))
                kf(params.uaux, params.qp.qe, params.B, params.fx, params.fy_t, params.metrics.Je,
                   params.ω, params.ω,
                   params.mesh.connijk, params.Minv,
                   params.mesh.ngl, params.mesh.ngl, params.neqs, lpert;
                   ndrange = (params.mesh.nelem * params.mesh.ngl, params.mesh.ngl),
                   workgroupsize = (params.mesh.ngl, params.mesh.ngl))
                KernelAbstractions.synchronize(backend)
                if (params.laguerre)
                    params.B_lag .= TFloat(0.0)
                    kf = filter_gpu_2d!(backend,(Int64(params.mesh.ngl), Int64(params.mesh.ngr)))
                    kf(params.uaux, params.qp.qe, params.B_lag, params.fx, params.fy_t_lag, params.metrics_lag.Je,
                       params.ω, params.ω_lag,
                       params.mesh.connijk_lag, params.Minv,
                       params.mesh.ngl, params.mesh.ngr, params.neqs, lpert;
                       ndrange = (params.mesh.nelem_semi_inf * params.mesh.ngl, params.mesh.ngr),
                       workgroupsize = (params.mesh.ngl, params.mesh.ngr))

                    KernelAbstractions.synchronize(backend)

                    params.B .+= params.B_lag
                end
                if (lpert)
                    params.uaux .= params.B
                else
                    params.uaux .= params.B .+ params.qp.qe
                end
                kf = uauxtou_gpu!(backend)
                kf(u,params.uaux,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
                   workgroupsize = (params.mesh.ngl,params.neqs))
                KernelAbstractions.synchronize(backend)
            end
            k = apply_boundary_conditions_gpu!(backend)
            k(@view(params.uaux[:,:]), @view(u[:]), params.qp.qe,
              params.mesh.x, params.mesh.y, TFloat(time),
              params.metrics.nx, params.metrics.ny,
              params.mesh.poin_in_bdy_edge,params.qbdy_gpu,
              params.mesh.ngl, TInt(params.neqs), params.mesh.npoin,lpert;
              ndrange = (params.mesh.nedges_bdy*params.mesh.ngl),
              workgroupsize = (params.mesh.ngl))

            KernelAbstractions.synchronize(backend)
            if (params.laguerre)

                k = apply_boundary_conditions_lag_gpu!(backend)
                k(@view(params.uaux[:,:]), @view(u[:]), params.qp.qe,
                  params.mesh.x, params.mesh.y, TFloat(time),
                  params.mesh.connijk_lag,
                  params.qbdy_lag_gpu,
                  params.mesh.ngl, params.mesh.ngr,
                  TInt(params.neqs), params.mesh.npoin, params.mesh.nelem_semi_inf,
                  params.inputs[:lperiodic_laguerre], lpert;
                  ndrange = (params.mesh.nelem_semi_inf*params.mesh.ngl,params.mesh.ngr),
                  workgroupsize = (params.mesh.ngl,params.mesh.ngr))

                KernelAbstractions.synchronize(backend)
            end

            k1(u,params.uaux,params.mesh.npoin,TInt(params.neqs);
               ndrange = (params.mesh.npoin,params.neqs),
               workgroupsize = (params.mesh.ngl,params.neqs))

            k = _build_rhs_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl)))

            k(params.RHS, params.uaux, params.qp.qe, params.mesh.x, params.mesh.y, params.mesh.connijk,
              params.metrics.dξdx, params.metrics.dξdy,
              params.metrics.dηdx, params.metrics.dηdy,
              params.metrics.Je,
              params.basis.dψ, params.ω, params.Minv, params.flux_gpu,
              params.source_gpu, params.mesh.ngl, TInt(params.neqs), PHYS_CONST,
              params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, lpert;
              ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl),
              workgroupsize = (params.mesh.ngl,params.mesh.ngl))

            KernelAbstractions.synchronize(backend)
            if (params.laguerre)
                params.RHS_lag .= TFloat(0.0)

                k_lag = _build_rhs_lag_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngr)))
                k_lag(params.RHS_lag, params.uaux, params.qp.qe,
                      params.mesh.x, params.mesh.y,
                      params.mesh.connijk_lag,
                      params.metrics_lag.dξdx, params.metrics_lag.dξdy,
                      params.metrics_lag.dηdx, params.metrics_lag.dηdy,
                      params.metrics_lag.Je,
                      params.basis.dψ, params.basis_lag.dψ, params.ω,
                      params.ω_lag, params.Minv, params.flux_lag_gpu, params.source_lag_gpu,
                      params.mesh.ngl, params.mesh.ngr, TInt(params.neqs), PHYS_CONST,
                      params.mesh.xmax, params.mesh.xmin, params.mesh.ymax, params.mesh.ymin, lpert;
                      ndrange = (params.mesh.nelem_semi_inf*params.mesh.ngl,params.mesh.ngr),
                      workgroupsize = (params.mesh.ngl,params.mesh.ngr))

                KernelAbstractions.synchronize(backend)

                @inbounds params.RHS .+= params.RHS_lag
                if (params.inputs[:lvisc])
                    params.RHS_visc_lag .= TFloat(0.0)
                    params.rhs_diffξ_el_lag .= TFloat(0.0)
                    params.rhs_diffη_el_lag .= TFloat(0.0)
                    params.source_lag_gpu .= TFloat(0.0)

                    k_diff_lag = _build_rhs_visc_lag_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngr)))
                    k_diff_lag(params.RHS_visc_lag,
                               params.rhs_diffξ_el_lag, params.rhs_diffη_el_lag,
                               params.uaux, params.qp.qe, params.source_lag_gpu,
                               params.mesh.x, params.mesh.y,
                               params.mesh.connijk_lag,
                               params.metrics_lag.dξdx, params.metrics_lag.dξdy,
                               params.metrics_lag.dηdx, params.metrics_lag.dηdy,
                               params.metrics_lag.Je, params.basis.dψ, params.basis_lag.dψ,
                               params.ω, params.ω_lag, params.Minv, params.visc_coeff,
                               params.mesh.ngl, params.mesh.ngr, TInt(params.neqs), PHYS_CONST, lpert;
                               ndrange = (params.mesh.nelem_semi_inf*params.mesh.ngl,params.mesh.ngr),
                               workgroupsize = (params.mesh.ngl,params.mesh.ngr))

                    @inbounds params.RHS .+= params.RHS_visc_lag

                end

            end

if (params.inputs[:lvisc])
    params.RHS_visc     .= TFloat(0.0)
    params.rhs_diffξ_el .= TFloat(0.0)
    params.rhs_diffη_el .= TFloat(0.0)
    params.source_gpu   .= TFloat(0.0)

    k = _build_rhs_diff_gpu_2D_v0!(backend, (Int64(params.mesh.ngl),Int64(params.mesh.ngl)))
    k(params.RHS_visc, params.rhs_diffξ_el, params.rhs_diffη_el,
      params.uaux, params.qp.qe, params.source_gpu,
      params.mesh.x, params.mesh.y, params.mesh.connijk,
      params.metrics.dξdx, params.metrics.dξdy,
      params.metrics.dηdx, params.metrics.dηdy,
      params.metrics.Je, params.basis.dψ, params.ω, params.Minv,
      params.visc_coeff, params.mesh.ngl, TInt(params.neqs), PHYS_CONST, lpert;
      ndrange = (params.mesh.nelem*params.mesh.ngl,params.mesh.ngl), workgroupsize = (params.mesh.ngl,params.mesh.ngl))
    KernelAbstractions.synchronize(backend)

    @inbounds params.RHS .+= params.RHS_visc
end
#@info maximum(params.RHS), maximum(params.RHS_lag), maximum(params.RHS_visc_lag)
DSS_global_RHS!(params.RHS, params.g_dss_cache, params.neqs)

k1 = RHStodu_gpu!(backend)
k1(params.RHS,du,params.mesh.npoin,TInt(params.neqs);ndrange = (params.mesh.npoin,params.neqs),
   workgroupsize = (params.mesh.ngl,params.neqs))

end
end
end # timer
end

function _build_rhs!(RHS, u, params, time)
@trixi_timeit timer() "build_rhs" begin
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

    inputs = params.inputs

    # comm = params.mesh.parts.comm
    comm = params.inputs.comm

    #-----------------------------------------------------------------------------------
    # Inviscid rhs:
    #-----------------------------------------------------------------------------------
    @trixi_timeit timer() " RESETRHSTOZERO "  resetRHSToZero_inviscid!(params)
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

    u2uaux!(params.uaux, u, params.neqs, params.mesh.npoin)

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
    
    #@code_warntype apply_boundary_conditions_dirichlet!(u, params.uaux, time, params.qp.qe,
    #@trixi_timeit timer() "apply DC boundary" apply_boundary_conditions_dirichlet!(u, params.uaux, time, params.qp.qe,
    apply_boundary_conditions_dirichlet!(u, params.uaux, time, params.qp.qe,
                                         params.mesh.coords,
                                         params.metrics.nx, params.metrics.ny, params.metrics.nz,
                                         params.mesh.npoin, params.mesh.npoin_linear,
                                         params.mesh.poin_in_bdy_edge, params.mesh.poin_in_bdy_face,
                                         params.mesh.nedges_bdy, params.mesh.nfaces_bdy, params.mesh.ngl,
                                         params.mesh.ngr, params.mesh.nelem_semi_inf, params.basis.ψ, params.basis.dψ,
                                         xmax, ymax, zmax, xmin, ymin, zmin, params.RHS, params.rhs_el, params.ubdy,
                                         params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem,
                                         params.mesh.bdy_edge_type, params.mesh.bdy_face_in_elem, params.mesh.bdy_face_type,
                                         params.mesh.connijk, params.metrics.Jef, params.S_face,
                                         params.S_flux, params.F_surf, params.M_surf_inv, params.M_edge_inv, params.Minv,
                                         params.mp.Tabs, params.mp.qn,
                                         params.ω, neqs, params.inputs, AD, SD)

    if (params.inputs[:lmoist])

        do_micro_physics!(params.mp.Tabs, params.mp.qn, params.mp.qc, params.mp.qi, params.mp.qr,
                          params.mp.qs, params.mp.qg, params.mp.Pr, params.mp.Ps, params.mp.Pg,
                          params.mp.S_micro, params.mp.qsatt, params.mesh.npoin,
                          params.uaux, @view(params.mesh.coords[:,end]),
                          params.qp.qe, SD, params.SOL_VARS_TYPE)

        if (params.inputs[:lprecip])
            compute_precipitation_derivatives!(params.mp.dqpdt, params.mp.dqtdt, params.mp.dhldt, params.mp.Pr, params.mp.Ps,
                                               params.mp.Pg, params.mp.Tabs, params.mp.qi,
                                               @view(params.uaux[:,1]), @view(params.qp.qe[:,1]),
                                               params.mesh.nelem, params.mesh.ngl, params.mesh.connijk, params.H,
                                               params.metrics, params.ω, params.basis.dψ, SD, params.SOL_VARS_TYPE)

            micro2rhs!(params.rhs_el,params.mp.dhldt, params.mp.dqtdt, params.mp.dqpdt, SD)
        end
        uaux2u!(u, params.uaux, params.neqs, params.mesh.npoin)
    end

    if(params.inputs[:lsaturation])
        saturation_adjustment(params.uaux, params.qp.qe, params.mesh.z, params.mesh.connijk,
                              params.mesh.nelem, params.mesh.ngl, neqs, params.thermo_params)

        uaux2u!(u, params.uaux, params.neqs, params.mesh.npoin)
    end



    @trixi_timeit timer() "inviscid_rhs_el!" inviscid_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, params.mesh.coords, lsource,
                     params.mp.S_micro, params.mp.qn, params.mp.flux_lw, params.mp.flux_sw, SD,
                     inputs.val_lsaturation)

    if inputs[:ladapt] == true
        DSS_nc_gather_rhs!(params.RHS, SD, QT, params.rhs_el,
                           params.mesh.non_conforming_facets,
                           params.mesh.non_conforming_facets_parents_ghost, params.cache_ghost_p,
                           params.q_el, params.q_el_pro, params.q_ghost_p,
                           params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPc_list_pg,
                           params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost,
                           params.mesh.pgip_local, ngl-1, neqs, params.interp)
    end
   @trixi_timeit timer() "DSS_rhs!" DSS_rhs!(params.RHS, params.rhs_el, params.mesh.connijk, nelem, ngl, neqs, SD, AD)

    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)

        resetRHSToZero_viscous!(params, SD)
        
        #Main.debug[] = (; u, params, connijk = params.mesh.connijk, qe = params.qp.qe, SD)
        #error()
        @trixi_timeit timer() "viscous_rhs_el!" viscous_rhs_el!(u, params, params.mesh.connijk, params.qp.qe, SD)
        
        if inputs[:ladapt] == true
            DSS_nc_gather_rhs!(params.RHS_visc, SD, QT, params.rhs_diff_el,
                               params.mesh.non_conforming_facets,
                               params.mesh.non_conforming_facets_parents_ghost, params.cache_ghost_p,
                               params.q_el, params.q_el_pro, params.q_ghost_p,
                               params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPc_list_pg,
                               params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost,
                               params.mesh.pgip_local, ngl-1, neqs, params.interp)
        end

        DSS_rhs!(params.RHS_visc, params.rhs_diff_el, params.mesh.connijk, nelem, ngl, neqs, SD, AD)
        params.RHS .+= params.RHS_visc
    end
    apply_boundary_conditions_neumann!(u, params.uaux, time, params.qp.qe,
                                       params.mesh.coords,
                                       params.metrics.nx, params.metrics.ny, params.metrics.nz,
                                       params.mesh.npoin, params.mesh.npoin_linear,
                                       params.mesh.poin_in_bdy_edge, params.mesh.poin_in_bdy_face,
                                       params.mesh.nedges_bdy, params.mesh.nfaces_bdy, params.mesh.ngl,
                                       params.mesh.ngr, params.mesh.nelem_semi_inf, params.basis.ψ, params.basis.dψ,
                                       xmax, ymax, zmax, xmin, ymin, zmin, params.RHS, params.rhs_el, params.ubdy,
                                       params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem,
                                       params.mesh.bdy_edge_type, params.mesh.bdy_face_in_elem, params.mesh.bdy_face_type,
                                       params.mesh.connijk, params.metrics.Jef, params.S_face,
                                       params.S_flux, params.F_surf, params.M_surf_inv, params.M_edge_inv, params.Minv,
                                       params.WM.τ_f, params.WM.wθ,
                                       params.mp.Tabs, params.mp.qn,
                                       params.ω, neqs, params.inputs, AD, SD)

    DSS_global_RHS!(params.RHS, params.g_dss_cache, params.neqs)

    #if (rem(time, Δt) == 0 && time > 0.0)
    if (time > 0.0)
        params.qp.qnm1 .= params.qp.qnm2
        params.qp.qnm2 .= params.uaux
    end

    Minv = params.Minv
    if Minv isa AbstractVector
        for ieq = 1:neqs
            divide_by_mass_matrix!(params.RHS, ieq, Minv, npoin, AD)
            if inputs[:ladapt] == true
                DSS_nc_scatter_rhs!(@view(params.RHS[:,ieq]), SD, QT,
                                    params.mesh.non_conforming_facets,
                                    params.mesh.non_conforming_facets_children_ghost, params.cache_ghost_c,
                                    params.q_el, params.q_el_pro, params.q_ghost_c,
                                    params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPp_list_cg,
                                    params.mesh.gip2ip, params.mesh.cgip_local, ngl-1, params.interp)
            end
        end
    else
        for ieq = 1:neqs
            divide_by_mass_matrix!(@view(params.RHS[:,ieq]), params.vaux, Minv, neqs, npoin, AD)
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
    end
end

function inviscid_rhs_el!(u, params,
                          connijk::Array{Int64,4},
                          qe::Matrix{Float64},
                          coords,
                          lsource, S_micro_vec, qn_vec,
                          flux_lw_vec, flux_sw_vec,
                          SD::NSD_1D, ::Val{false})

    u2uaux!(params.uaux, u, params.neqs, params.mesh.npoin)

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
                             @view(params.uprimitive[:,i]),
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
                             params.uprimitive,
                             params.F, params.S,
                             params.rhs_el,
                             iel, params.CL, params.QT, SD, params.AD)

    end
end


@inline function _expansion_inviscid_KEP!(u, neqs, ngl,
                                  dψ, ω,
                                  F, S, D,
                                  rhs_el, uilgl, iel,
                                  ::CL, QT::Inexact,
                                  SD::NSD_1D, AD::ContGal,
                                  uaux, connijk, el,
                                  volume_flux_type)

    for i = 1:ngl
        ip = connijk[el,i,1]
        du_i = zeros(neqs)

        for j = 1:ngl
            jp = connijk[el,j,1]
            f_ij = user_volume_flux(uaux[ip,:], uaux[jp,:], volume_flux_type)
            for ieq = 1:neqs
                du_i[ieq] += 2.0 *  dψ[j, i] * f_ij[ieq]
            end
        end

        for ieq = 1:neqs
            rhs_el[iel, i, ieq] -=  ω[i] *du_i[ieq] - ω[i] * S[i, ieq]
        end
    end
end


@inline function _expansion_inviscid_KEP!(u, neqs, ngl, dψ, ω,
                                  F, G, S,
                                  Je,
                                  dξdx, dξdy,
                                  dηdx, dηdy,
                                  rhs_el, iel,
                                  ::CL, QT::Inexact,
                                  SD::NSD_2D, AD::ContGal,
                                  uaux, fluxaux,
  				     dFdx, dFdxi, dFdeta,
				     dGdy, dGdxi, dGdeta,
				  connijk, volume_flux_type)

    for j=1:ngl
        for i=1:ngl
            ip = connijk[iel,i,j]
            ωJac = ω[i]*ω[j]*Je[i, j, iel]
            @. dFdxi = 0
	    @. dFdeta = 0
	    @. dGdxi = 0
	    @. dGdeta = 0
            for k = 1:ngl
                kjp = connijk[iel,k, j]
        	ikp = connijk[iel,i, k]

	F_ik, G_ik = flux_turbo(@view(fluxaux[ip,:]), @view(fluxaux[ikp,:]), volume_flux_type)
	F_kj, G_kj = flux_turbo(@view(fluxaux[ip,:]), @view(fluxaux[kjp,:]), volume_flux_type)
                 @. dFdxi += 2 * dψ[k,i]*F_kj
                 @. dFdeta += 2 * dψ[k,j]*F_ik
                 @. dGdxi += 2 * dψ[k,i]*G_kj
                 @. dGdeta += 2 * dψ[k,j]*G_ik
            end
            dξdx_ij = dξdx[i, j, iel]
            dξdy_ij = dξdy[i, j, iel]
            dηdx_ij = dηdx[i, j, iel]
            dηdy_ij = dηdy[i, j, iel]

             @. dFdx = dFdxi*dξdx_ij + dFdeta*dηdx_ij
  	     @. dGdy = dGdxi*dξdy_ij + dGdeta*dηdy_ij
            for ieq=1:neqs
                rhs_el[iel,i,j,ieq] -=  ωJac*((dFdx[ieq] + dGdy[ieq]) - S[i,j,ieq])
            end
        end
    end
end

function inviscid_rhs_el!(u, params,
                          connijk::Array{Int64,4},
                          qe::Matrix{Float64},
                          coords,
                          lsource, S_micro_vec, qn_vec, flux_lw_vec,
                          flux_sw_vec, SD::NSD_2D, ::Val)

    ngl   = params.mesh.ngl
    nelem = params.mesh.nelem

    xmin = params.xmin; xmax = params.xmax; ymax = params.ymax

    lkep = params.inputs[:lkep]

    for iel = 1:nelem

        for j = 1:ngl, i=1:ngl
            ip = connijk[iel,i,j]

            user_primitives!(@view(params.uaux[ip,:]),@view(qe[ip,:]),@view(params.uprimitive[:,i,j]), params.SOL_VARS_TYPE)
            if lkep
         user_fluxaux!(@view(params.fluxaux[ip,:]),
                              SD,
                              @view(params.uaux[ip,:]),
                              params.SOL_VARS_TYPE,
                              params.volume_flux)
            else
                user_flux!(@view(params.F[i,j,:]), @view(params.G[i,j,:]), SD,
                           @view(params.uaux[ip,:]),
                           @view(qe[ip,:]),
                           params.mesh,
                           params.CL, params.SOL_VARS_TYPE;
                           neqs=params.neqs, ip=ip)
            end

            if lsource
          @trixi_timeit timer() "user_source!"      user_source!(@view(params.S[i,j,:]),
                             @view(params.uaux[ip,:]),
                             @view(qe[ip,:]),
                             params.mesh.npoin, params.CL, params.SOL_VARS_TYPE;
                             neqs=params.neqs,
                             x=coords[ip,1], y=coords[ip,2], ymax=ymax)

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

        if lkep
           _expansion_inviscid_KEP!(u,
                                     params.neqs, params.mesh.ngl,
                                     params.basis.dψ, params.ω,
                                     params.F, params.G, params.S,
                                     params.metrics.Je,
                                     params.metrics.dξdx, params.metrics.dξdy,
                                     params.metrics.dηdx, params.metrics.dηdy,
                                     params.rhs_el, iel, params.CL, params.QT, SD,
                                     params.AD, params.uaux, params.fluxaux,
				     params.dFdx, params.dFdxi, params.dFdeta,
				     params.dGdy, params.dGdxi, params.dGdeta,
				     connijk,
                                     params.volume_flux)
        else
            _expansion_inviscid!(u,
                                 params.neqs, params.mesh.ngl,
                                 params.basis.dψ, params.ω,
                                 params.uprimitive,
                                 params.F, params.G, params.S,
                                 params.metrics.Je,
                                 params.metrics.dξdx, params.metrics.dξdy,
                                 params.metrics.dηdx, params.metrics.dηdy,
                                 params.rhs_el, iel, params.CL, params.QT, SD, params.AD)
        end


    end
end

function inviscid_rhs_el!(u, params,
                          connijk::Array{Int64,4},
                          qe::Matrix{Float64},
                          coords,
                          lsource, S_micro_vec, qn_vec,
                          flux_lw_vec, flux_sw_vec,
                          SD::NSD_3D, ::Val{false})

    nelem = params.mesh.nelem
    ngl   = params.mesh.ngl

    u2uaux!(params.uaux, u, params.neqs, params.mesh.npoin)
    xmin = params.xmin; xmax = params.xmax; zmax = params.zmax
    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i=1:ngl

            ip = connijk[iel,i,j,k]

            user_primitives!(@view(params.uaux[ip,:]),@view(qe[ip,:]),@view(params.uprimitive[:,i,j,k]), params.SOL_VARS_TYPE)


            user_flux!(@view(params.F[i,j,k,:]),
                       @view(params.G[i,j,k,:]),
                       @view(params.H[i,j,k,:]),
                       @view(params.uaux[ip,:]),
                       @view(qe[ip,:]),
                       params.mesh,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip)

            if lsource
                user_source!(@view(params.S[i,j,k,:]),
                             @view(params.uaux[ip,:]),
                             @view(qe[ip,:]),
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
                                              S_micro, qn, flux_lw, flux_sw, PHYS_CONST,
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
                             params.uprimitive,
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





function inviscid_rhs_el!(u, params, connijk::Array{Int64,4}, qe::Matrix{Float64}, coords, lsource, S_micro_vec, qn_vec, flux_lw_vec, flux_sw_vec, SD::NSD_3D, ::Val{true})

    nelem = params.mesh.nelem
    ngl   = params.mesh.ngl

    u2uaux!(params.uaux, u, params.neqs, params.mesh.npoin)
    xmin = params.xmin; xmax = params.xmax; zmax = params.zmax

    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i=1:ngl

            ip = connijk[iel,i,j,k]

            user_flux!(@view(params.F[i,j,k,:]),
                       @view(params.G[i,j,k,:]),
                       @view(params.H[i,j,k,:]),
                       @view(params.uaux[ip,:]),
                       @view(qe[ip,:]),
                       params.mesh, params.thermo_params,
                       params.CL, params.SOL_VARS_TYPE;
                       neqs=params.neqs, ip=ip,
                       x=coords[ip,1], y=coords[ip,2], z=coords[ip,3])

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
                                              S_micro, qn, flux_lw, flux_sw, PHYS_CONST,
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



function viscous_rhs_el!(u, params, connijk::Array{Int64,4}, qe::Matrix{Float64}, SD::NSD_1D)

    Δ = params.mesh.Δeffective_l

    nelem = params.mesh.nelem
    ngl   = params.mesh.ngl
    neqs  = params.neqs

    for iel=1:nelem

        for i=1:ngl
            ip = connijk[iel,i]
           user_primitives!(@view(params.uaux[ip,:]), @view(qe[ip,:]), @view(params.uprimitive[:,i]), params.SOL_VARS_TYPE)
        end

        for ieq = 1:neqs
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

@inbounds function viscous_rhs_el!(u, params, connijk, qe, SD::NSD_2D)

    @unpack mesh, neqs, inputs, uaux, uprimitive, SOL_VARS_TYPE, rhs_diffξ_el, rhs_diffη_el, visc_coeff, ω, metrics, rhs_el, ω, QT, VT, AD, basis, gradient_dxi, gradient_deta, gradient_dx, gradient_dy, dx_flux, dy_flux, rhs_diff_el = params

    @unpack entropy_variables = inputs
    @unpack Je, dξdx, dξdy, dηdx, dηdy = metrics
    @unpack dψ = basis
    @unpack nelem, ngl = mesh
    Δ = mesh.Δeffective_l

    # Add call to a function that takes all the above arguments as inputs
    viscous_rhs_el_type_stable!(u, connijk, neqs, nelem, ngl, dψ, ω,
                                  rhs_diffξ_el, rhs_diffη_el,
                                  uaux, qe,
                                  uprimitive,
                                  SOL_VARS_TYPE,
                                  visc_coeff,
                                  Je,
                                  dξdx, dξdy,
                                  dηdx, dηdy,
                                  inputs,
                                  gradient_dxi, gradient_deta,
                                  gradient_dx, gradient_dy, dx_flux, dy_flux,
                                  rhs_el, QT, VT, SD, AD, Δ, rhs_diff_el, entropy_variables)
end

@inbounds function viscous_rhs_el_type_stable!(u, connijk, neqs, nelem, ngl, dψ, ω,
                                  rhs_diffξ_el, rhs_diffη_el,
                                  uaux, qe,
                                  uprimitive,
                                  SOL_VARS_TYPE,
                                  visc_coeff,
                                  Je,
                                  dξdx, dξdy,
                                  dηdx, dηdy,
                                  inputs,
                                  gradient_dxi, gradient_deta,
                                  gradient_dx, gradient_dy, dx_flux, dy_flux,
                                  rhs_el, QT, VT, SD, AD, Δ, rhs_diff_el, entropy_variables)


    if entropy_variables
        # Compute the u_transformed everywhere and store in uprimitive
        for iel=1:nelem
            for j = 1:ngl, i=1:ngl
                ip = connijk[iel,i,j]
                uaux_node = get_node_vars_4(uaux, ip)
                # qe_node = get_node_vars_4(qe, ip)
                # uprimitive_node = get_node_vars_4(uprimitive, i, j)
                @views user_primitives!(uaux[ip,:], qe[ip,:], uprimitive[:,i,j], SOL_VARS_TYPE)
            end
            _expansion_visc_navierstokes!(rhs_diffξ_el,
                            rhs_diffη_el,
                            uprimitive,
                            visc_coeff,
                            ω,
                            ngl,
                            dψ,
                            Je,
                            dξdx, dξdy,
                            dηdx, dηdy,
                            inputs, rhs_el,
                            iel, neqs, gradient_dxi, gradient_deta,
                            gradient_dx, gradient_dy, dx_flux, dy_flux,
                            QT, VT, SD, AD; Δ=Δ)
        end


    else
        for iel=1:nelem
            
            @inbounds for j = 1:ngl, i=1:ngl
                ip = connijk[iel,i,j]
                # uaux_node = @inline get_node_vars_4(uaux, ip)
                # qe_node = @inline get_node_vars_4(qe, ip)
                # uprimitive_node = @inline get_node_vars_4(uprimitive, i, j)
                # user_primitives_node = user_primitives(uaux_node, qe_node, uprimitive_node, SOL_VARS_TYPE)
                # set_node_vars_4!(uaux, user_primitives_node, ip)
                @views user_primitives!(uaux[ip,:],qe[ip,:],uprimitive[:,i,j], SOL_VARS_TYPE)
            end

            for ieq = 1:neqs
                _expansion_visc!(rhs_diffξ_el,
                                rhs_diffη_el,
                                uprimitive,
                                visc_coeff,
                                ω,
                                ngl,
                                dψ,
                                Je,
                                dξdx, dξdy,
                                dηdx, dηdy,
                                inputs, rhs_el,
                                iel, ieq,
                                QT, VT, SD, AD; Δ=Δ)
            end

        end
    end

    @. rhs_diff_el = rhs_diffξ_el + rhs_diffη_el
end


function viscous_rhs_el!(u, params, connijk::Array{Int64,4}, qe::Matrix{Float64}, SD::NSD_3D)

    Δ = params.mesh.Δeffective_l

    nelem = params.mesh.nelem
    ngl   = params.mesh.ngl
    neqs  = params.neqs

    fill!(params.μ_max,    zero(params.T))

    for iel=1:nelem

        for k = 1:ngl, j = 1:ngl, i=1:ngl
            ip = connijk[iel,i,j,k]

            user_primitives!(@view(params.uaux[ip,:]),
                             @view(qe[ip,:]),
                             @view(params.uprimitive[:,i,j,k]),
                             params.SOL_VARS_TYPE)
        end


        for ieq = 1:neqs
            _expansion_visc!(params.rhs_diffξ_el,
                             params.rhs_diffη_el,
                             params.rhs_diffζ_el,
                             params.uprimitive,
                             params.visc_coeff,
                             params.ω,
                             params.mp.Tabs,
                             params.mp.qn,
                             params.mp.qsatt,
                             params.uaux,
                             params.mesh.ngl,
                             params.basis.dψ,
                             params.metrics.Je,
                             params.metrics.dξdx, params.metrics.dξdy, params.metrics.dξdz,
                             params.metrics.dηdx, params.metrics.dηdy, params.metrics.dηdz,
                             params.metrics.dζdx, params.metrics.dζdy, params.metrics.dζdz,
                             params.inputs, params.rhs_el, iel, ieq, params.mesh.connijk,
                             params.mesh.coords,
                             params.mesh.poin_in_bdy_face, params.mesh.elem_to_face,
                             params.mesh.bdy_face_type,
                             params.μ_max,
                             params.QT, params.VT, SD, params.AD; Δ=Δ)

        end
    end

    params.rhs_diff_el .= @views (params.rhs_diffξ_el .+ params.rhs_diffη_el .+ params.rhs_diffζ_el)

end


function _expansion_inviscid!(u, params, iel, ::CL, QT::Inexact, SD::NSD_1D, AD::FD)

    ngl   = params.mesh.ngl
    neqs  = params.neqs
    npoin = params.mesh.npoin

    for ieq = 1:neqs
        for i = 1:ngl
            ip = params.mesh.connijk[iel,i,1]
            if (ip < npoin)
                params.RHS[ip,ieq] = 0.5*(u[ip+1] - u[ip])/(params.mesh.Δx[ip])
            end
        end
    end
    nothing
end


function _expansion_inviscid!(u, neqs, ngl,
                              dψ, ω,
                              uprimitive,
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

function _expansion_inviscid!(u, neqs, ngl,
                              dψ, ω,
                              uprimitive,
                              F, G, S,
                              Je,
                              dξdx, dξdy,
                              dηdx, dηdy,
                              rhs_el, iel,
                              ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for ieq=1:neqs
        for j=1:ngl
            ωj = ω[j]
            for i=1:ngl

                @inbounds begin
                    Jeij = Je[i, j, iel]
                    ωJac = ω[i]*ωj*Jeij

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
                    dξdx_ij = dξdx[i, j, iel]
                    dξdy_ij = dξdy[i, j, iel]
                    dηdx_ij = dηdx[i, j, iel]
                    dηdy_ij = dηdy[i, j, iel]

                    dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
                    dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

                    rhs_el[iel,i,j,ieq] -=  ωJac*((dFdx + dGdy) - S[i,j,ieq])
                end
            end
        end
    end
end

function _expansion_inviscid!(u, neqs, ngl,
                              dψ, ω,
                              uprimitive,
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

                ωj = ω[j]
                ωk = ω[k]
                ωjk = ωj * ωk

                for i=1:ngl

                    @inbounds begin
                        Je_ijk = Je[i, j, k, iel]
                        ωJac = ω[i] * ωjk * Je_ijk

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
                        dξdx_ij = dξdx[i, j, k, iel]
                        dξdy_ij = dξdy[i, j, k, iel]
                        dξdz_ij = dξdz[i, j, k, iel]

                        dηdx_ij = dηdx[i, j, k, iel]
                        dηdy_ij = dηdy[i, j, k, iel]
                        dηdz_ij = dηdz[i, j, k, iel]

                        dζdx_ij = dζdx[i, j, k, iel]
                        dζdy_ij = dζdy[i, j, k, iel]
                        dζdz_ij = dζdz[i, j, k, iel]

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
end

function _expansion_inviscid!(u, params, iel, ::CL, QT::Exact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(u, neqs, ngl,
                              dψ, ω,
                              uprimitive,
                              F, G, S,
                              Je,
                              dξdx, dξdy,
                              dηdx, dηdy,
                              rhs_el, iel,
                              ::CL, QT::Exact, SD::NSD_2D, AD::ContGal)

    N    = ngl
    Q    = N + 1

    for ieq=1:neqs
        for l=1:Q
            ωl = ω[l]
            for k=1:Q
                @inbounds begin
                    Je_kl = Je[k, l, iel]
                    ωJac = ω[k] * ωl * Je_kl

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

                    dξdx_kl = params.metrics.dξdx[k, l, iel]
                    dξdy_kl = params.metrics.dξdy[k, l, iel]
                    dηdx_kl = params.metrics.dηdx[k, l, iel]
                    dηdy_kl = params.metrics.dηdy[k, l, iel]
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
end

function _expansion_inviscid!(u, params, iel, ::NCL, QT::Inexact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(u, neqs, ngl,
                              dψ, ω,
                              uprimitive,
                              F, G, S,
                              Je,
                              dξdx, dξdy,
                              dηdx, dηdy,
                              rhs_el, iel,
                              ::NCL, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for ieq=1:neqs
        for j=1:ngl
            ωj = ω[j]
            for i=1:ngl

                @inbounds begin
                    Je_ij = Je[i, j, iel]
                    ωJac  = ω[i]*ωj*Je_ij

                    dFdξ = 0.0; dFdη = 0.0
                    dGdξ = 0.0; dGdη = 0.0
                    dpdξ = 0.0; dpdη = 0.0
                    for k = 1:ngl
                        dFdξ += dψ[k,i]*F[k,j,ieq]
                        dFdη += dψ[k,j]*F[i,k,ieq]

                        dGdξ += dψ[k,i]*G[k,j,ieq]
                        dGdη += dψ[k,j]*G[i,k,ieq]

                        dpdξ += dψ[k,i]*uprimitive[neqs+1,k,j]
                        dpdη += dψ[k,j]*uprimitive[neqs+1,i,k]
                    end
                    dξdx_ij = dξdx[i, j, iel]
                    dξdy_ij = dξdy[i, j, iel]
                    dηdx_ij = dηdx[i, j, iel]
                    dηdy_ij = dηdy[i, j, iel]

                    dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
                    dFdy = dFdξ*dξdy_ij + dFdη*dηdy_ij

                    dGdx = dGdξ*dξdx_ij + dGdη*dηdx_ij
                    dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

                    dpdx = dpdξ*dξdx_ij + dpdη*dηdx_ij
                    dpdy = dpdξ*dξdy_ij + dpdη*dηdy_ij

                    ρij = uprimitive[1,i,j]
                    uij = uprimitive[2,i,j]
                    vij = uprimitive[3,i,j]

                    if (ieq == 1)
                        auxi = ωJac*(dFdx + dGdy)
                    elseif(ieq == 2)
                        auxi = ωJac*(uij*dFdx + vij*dGdy + dpdx/ρij)
                    elseif(ieq == 3)
                        auxi = ωJac*(uij*dFdx + vij*dGdy + dpdy/ρij - S[i,j,ieq])
                    elseif(ieq == 4)
                        auxi = ωJac*(uij*dFdx + vij*dGdy)
                    end

                    rhs_el[iel,i,j,ieq] -= auxi
                end
            end
        end
    end
end


function _expansion_inviscid!(u, params, iel, ::NCL, QT::Exact, SD::NSD_2D, AD::FD) nothing end

function _expansion_inviscid!(u, params, iel, ::NCL, QT::Exact, SD::NSD_2D, AD::ContGal)

    N = params.mesh.ngl
    Q = N + 1

    for l=1:Q
        ωl = ω[l]
        for k=1:Q

            @inbounds begin
                ωJac = ω[k]*ωl*Je[k, l, iel]

                dρudξ = 0.0; dρudη = 0.0
                dρvdξ = 0.0; dρvdη = 0.0
                dudξ  = 0.0; dudη  = 0.0
                dvdξ  = 0.0; dvdη  = 0.0
                dθdξ  = 0.0; dθdη  = 0.0
                dpdξ  = 0.0; dpdη  = 0.0

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

                        dudξ += dψmk_ψnl*params.uprimitive[2,m,n]
                        dudη +=  ψmk_dψnl*params.uprimitive[2,m,n]

                        dvdξ += dψmk_ψnl*params.uprimitive[3,m,n]
                        dvdη +=  ψmk_dψnl*params.uprimitive[3,m,n]

                        dθdξ += dψmk_ψnl*params.uprimitive[4,m,n]
                        dθdη +=  ψmk_dψnl*params.uprimitive[4,m,n]

                        dpdξ += dψmk_ψnl*params.uprimitive[params.neqs+1,m,n]
                        dpdη +=  ψmk_dψnl*params.uprimitive[params.neqs+1,m,n]

                        ρkl += ψmk*ψnl*params.uprimitive[1,m,n]
                        ukl += ψmk*ψnl*params.uprimitive[2,m,n]
                        vkl += ψmk*ψnl*params.uprimitive[3,m,n]
                        Skl += ψmk*ψnl*params.S[m,n,3]
                    end
                end

                dξdx_kl = params.metrics.dξdx[k, l, iel]
                dξdy_kl = params.metrics.dξdy[k, l, iel]
                dηdx_kl = params.metrics.dηdx[k, l, iel]
                dηdy_kl = params.metrics.dηdy[k, l, iel]

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
end


@inline function _expansion_visc!(rhs_diffξ_el, uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je, dξdx, inputs, rhs_el, iel, ieq,
                          QT::Inexact, VT::AV, SD::NSD_1D, ::ContGal; Δ=1.0)

    for k = 1:ngl
        ωJac = ω[k]*Je[k, iel]

        dqdξ = 0.0
        @turbo for ii = 1:ngl
            dqdξ += dψ[ii,k]*uprimitiveieq[ieq,ii]
        end

        dξdx_kl = dqdξ*dξdx[k, iel]
        dqdx = visc_coeffieq[ieq]*dξdx_kl

        ∇ξ∇u_kl = dξdx[k, iel]*dqdx*ωJac

        @turbo for i = 1:ngl
            dhdξ_ik = dψ[i,k]

            rhs_diffξ_el[iel,i,ieq] -= dhdξ_ik * ∇ξ∇u_kl
        end
    end
end


@inline function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeffieq, ω,
                          mesh, basis, metrics, inputs, rhs_el, iel, ieq,
                          QT::Inexact, VT, SD::NSD_2D, ::FD)
    nothing
end

@inline function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, ieq,
                          QT::Inexact, VT::AV, SD::NSD_2D, ::ContGal; Δ=1.0)

    for l = 1:ngl
        ωl = ω[l]
        for k = 1:ngl

            @inbounds begin
                Jekl = Je[k, l, iel]
                ωJac = ω[k]*ωl*Jekl

                dqdξ = 0.0
                dqdη = 0.0
                @turbo for ii = 1:ngl
                    dqdξ += dψ[ii,k]*uprimitiveieq[ieq,ii,l]
                    dqdη += dψ[ii,l]*uprimitiveieq[ieq,k,ii]
                end
                dξdx_kl = dξdx[k, l, iel]
                dξdy_kl = dξdy[k, l, iel]
                dηdx_kl = dηdx[k, l, iel]
                dηdy_kl = dηdy[k, l, iel]

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
end

# viscous RHS 2D
# SMAG FUNCTION
@inline function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, ieq,
                          QT::Inexact, VT, SD::NSD_2D, ::ContGal; Δ=1.0, vargs...)

    Sc_t      = PHYS_CONST.Sc_t
    Δ2        = Δ^2

    # Determine if this is a momentum equation
    is_u_momentum  = (ieq == 2)
    is_v_momentum  = (ieq == 3)
    is_temperature = (ieq == 4)

    for l = 1:ngl
        ωl = ω[l]
        for k = 1:ngl

            @inbounds begin
                Je_kl = Je[k, l, iel]
                ωJac  = ω[k]*ωl*Je_kl

                # Quantities for Smagorinsky
                dudξ = 0.0; dudη = 0.0
                dvdξ = 0.0; dvdη = 0.0
		## Computing the gradients
                @turbo for ii = 1:ngl
                    dudξ += dψ[ii,k]*uprimitiveieq[2,ii,l]
                    dudη += dψ[ii,l]*uprimitiveieq[2,k,ii]
                    dvdξ += dψ[ii,k]*uprimitiveieq[3,ii,l]
                    dvdη += dψ[ii,l]*uprimitiveieq[3,k,ii]
                end
                dξdx_kl = dξdx[k, l, iel]
                dξdy_kl = dξdy[k, l, iel]
                dηdx_kl = dηdx[k, l, iel]
                dηdy_kl = dηdy[k, l, iel]

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
                                                         uprimitiveieq[1,k,l],
                                                         dudx, dvdy, dudy, dvdx,
                                                         PHYS_CONST, Δ2,
                                                         inputs,
                                                         VT, SD)

                    τ_xx = 2.0 * effective_viscosity * dudx - (2.0/3.0) * effective_viscosity * div_u
                    τ_xy = effective_viscosity * (dudy + dvdx)
                    flux_x = τ_xx
                    flux_y = τ_xy


                elseif is_v_momentum
                    # USE EFFECTIVE VISCOSITY
                    effective_viscosity =  SGS_diffusion(visc_coeffieq, ieq,
                                                         uprimitiveieq[1,k,l],
                                                         dudx, dvdy, dudy, dvdx,
                                                         PHYS_CONST, Δ2,
                                                         inputs,
                                                         VT, SD)

                    τ_xy = effective_viscosity * (dudy + dvdx)
                    τ_yy = 2.0 * effective_viscosity * dvdy - (2.0/3.0) * effective_viscosity * div_u
                    flux_x = τ_xy
                    flux_y = τ_yy

                elseif is_temperature
                    # USE EFFECTIVE DIFFUSIVITY
                    effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                          uprimitiveieq[1,k,l],
                                                          dudx, dvdy, dudy, dvdx,
                                                          PHYS_CONST, Δ2,
                                                          inputs,
                                                          VT, SD)

                    # Compute temperature gradient
                    dθdξ = 0.0; dθdη = 0.0
                    @turbo for ii = 1:ngl
                        dθdξ += dψ[ii,k]*uprimitiveieq[ieq,ii,l]
                        dθdη += dψ[ii,l]*uprimitiveieq[ieq,k,ii]
                    end

                    dθdx = dθdξ*dξdx_kl + dθdη*dηdx_kl
                    dθdy = dθdξ*dξdy_kl + dθdη*dηdy_kl

                    flux_x = effective_diffusivity * dθdx
                    flux_y = effective_diffusivity * dθdy

                else
                    # Other scalars (use appropriate Schmidt number)
                    # USE EFFECTIVE DIFFUSIVITY
                    effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                          uprimitiveieq[1,k,l],
                                                          dudx, dvdy, dudy, dvdx,
                                                          PHYS_CONST, Δ2,
                                                          inputs,
                                                          VT, SD)

                    # Compute temperature gradient
                    dqdξ = 0.0; dqdη = 0.0
                    @turbo for ii = 1:ngl
                        dqdξ += dψ[ii,k]*uprimitiveieq[ieq,ii,l]
                        dqdη += dψ[ii,l]*uprimitiveieq[ieq,k,ii]
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
end

@inline @inbounds function _expansion_visc_navierstokes!(rhs_diffξ_el, rhs_diffη_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          ngl, dψ, Je,
                          dξdx, dξdy,
                          dηdx, dηdy,
                          inputs, rhs_el,
                          iel, neqs,
			  gradient_dxi, gradient_deta,
			  gradient_dx, gradient_dy, dx_flux, dy_flux,
                          QT::Inexact, VT, SD::NSD_2D, ::ContGal; Δ=1.0, vargs...)

    Sc_t      = PHYS_CONST.Sc_t
    Δ2        = Δ^2
    
    for l = 1:ngl
        ωl = ω[l]
        for k = 1:ngl

            @inbounds begin
                Je_kl = Je[k, l, iel]
                ωJac  = ω[k]*ωl*Je_kl
                @. gradient_dxi = 0
                @. gradient_deta = 0
		## Computing the gradients
		for var in  1:neqs
                for ii = 1:ngl
		    gradient_dxi[var] += dψ[ii,k]*uprimitiveieq[var,ii,l]
		    gradient_deta[var] += dψ[ii,l]*uprimitiveieq[var,k,ii]
                end
		end
                dξdx_kl = dξdx[k, l, iel]
                dξdy_kl = dξdy[k, l, iel]
                dηdx_kl = dηdx[k, l, iel]
                dηdy_kl = dηdy[k, l, iel]

                @. gradient_dx = gradient_dxi*dξdx_kl + gradient_deta*dηdx_kl
                @. gradient_dy = gradient_dxi*dξdy_kl + gradient_deta*dηdy_kl

		## TODO: Compute parabolic fluxes
		@views flux_x = flux_parabolic(uprimitiveieq[:,k,l], (gradient_dx, gradient_dy), 1, visc_coeffieq, inputs, Δ2)
		@views flux_y = flux_parabolic(uprimitiveieq[:,k,l], (gradient_dx, gradient_dy), 2, visc_coeffieq, inputs,Δ2)
		## FIX: reference or physical and arrays
                @. dx_flux = (dξdx_kl*flux_x + dξdy_kl*flux_y)*ωJac
                @. dy_flux = (dηdx_kl*flux_x + dηdy_kl*flux_y)*ωJac
                @turbo for i = 1:ngl
                    dhdξ_ik = dψ[i,k]
                    dhdη_il = dψ[i,l]
                   for ieq in 1:neqs
		    rhs_diffξ_el[iel,i,l,ieq] -= dhdξ_ik * dx_flux[ieq]
		    rhs_diffη_el[iel,k,i,ieq] -= dhdη_il * dy_flux[ieq]
	    	  end
                end
            end
        end
    end
end

@inline function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          Tabs, qn, qs,
                          uaux,
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
                          μsgs,
                          QT::Inexact, VT::AV, SD::NSD_3D, ::ContGal; Δ=1.0)
    conn_el = @view connijk[iel,:,:,:]
    lsponge = inputs[:lsponge]
    zs      = inputs[:zsponge]
    for m = 1:ngl
        for l = 1:ngl

            ωl = ω[l]
            ωm = ω[m]
            ωlm = ωl * ωm

            for k = 1:ngl

                @inbounds begin
                    Je_klm = Je[k, l, m, iel]
                    ωJac   = ω[k] * ωlm * Je_klm
                    ip     = conn_el[k,l,m]
                    z      = coords[ip,3]

                    σμ     = 1.0
                    if (z > zs) && (ieq > 4)
                        Z = (z - zs) / (25000. - zs)
                        # Formula: 1 - (10*X^3 - 15*X^4 + 6*X^5)
                        σμ = 1 - (Z^3 * (10.0 + Z * (-15.0 + Z * 6.0)))
                    end

                    dqdξ = 0.0
                    dqdη = 0.0
                    dqdζ = 0.0
                    @turbo for ii = 1:ngl
                        dqdξ += dψ[ii,k]*uprimitiveieq[ieq,ii,l,m]
                        dqdη += dψ[ii,l]*uprimitiveieq[ieq,k,ii,m]
                        dqdζ += dψ[ii,m]*uprimitiveieq[ieq,k,l,ii]
                    end
                    dξdx_klm = dξdx[k, l, m, iel]
                    dξdy_klm = dξdy[k, l, m, iel]
                    dξdz_klm = dξdz[k, l, m, iel]

                    dηdx_klm = dηdx[k, l, m, iel]
                    dηdy_klm = dηdy[k, l, m, iel]
                    dηdz_klm = dηdz[k, l, m, iel]

                    dζdx_klm = dζdx[k, l, m, iel]
                    dζdy_klm = dζdy[k, l, m, iel]
                    dζdz_klm = dζdz[k, l, m, iel]

                    auxi = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
                    dqdx = visc_coeffieq[ieq]*auxi

                    auxi = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
                    dqdy = visc_coeffieq[ieq]*auxi

                    auxi = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
                    dqdz = visc_coeffieq[ieq]*auxi

                    ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz)*ωJac * σμ
                    ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz)*ωJac * σμ
                    ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz)*ωJac * σμ

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
end



@inline function _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                          uprimitiveieq, visc_coeffieq, ω,
                          Tabs, qn, qs,
                          uaux,
                          ngl, dψ, Je,
                          dξdx, dξdy, dξdz,
                          dηdx, dηdy, dηdz,
                          dζdx, dζdy, dζdz,
                          inputs, rhs_el,
                          iel, ieq, connijk,
                          coords,
                          poin_in_bdy_face, elem_to_face, bdy_face_type,
                          μ_max,
                          QT::Inexact, VT, SD::NSD_3D, ::ContGal; Δ=1.0)

    Δ2 = Δ^2

    # Determine equation type (indices shifted for 3D)
    is_u_momentum  = (ieq == 2)
    is_v_momentum  = (ieq == 3)
    is_w_momentum  = (ieq == 4)
    is_temperature = (ieq == 5)
    conn_el        = @view connijk[iel,:,:,:]
    μ_max_ieq      = μ_max[ieq]

    lsponge = inputs[:lsponge]
    zs      = inputs[:zsponge]

    for m = 1:ngl
        for l = 1:ngl

            ωl = ω[l]
            ωm = ω[m]
            ωlm = ωl * ωm

            for k = 1:ngl

                ip     = conn_el[k,l,m]
                z      = coords[ip,3]

                σμ     = 1.0
                if (z > zs) && (ieq > 4)
                    Z = (z - zs) / (25000. - zs)
                    # Formula: 1 - (10*X^3 - 15*X^4 + 6*X^5)
                    σμ = 1 - (Z^3 * (10.0 + Z * (-15.0 + Z * 6.0)))
                end
                @inbounds begin
                    Je_klm = Je[k, l, m, iel]
                    ωJac = ω[k] * ωlm * Je_klm

                    # ===== Compute all velocity gradients =====
                    # u-velocity gradients
                    dudξ = 0.0; dudη = 0.0; dudζ = 0.0
                    @turbo for ii = 1:ngl
                        dudξ += dψ[ii,k]*uprimitiveieq[2,ii,l,m]
                        dudη += dψ[ii,l]*uprimitiveieq[2,k,ii,m]
                        dudζ += dψ[ii,m]*uprimitiveieq[2,k,l,ii]
                    end

                    # v-velocity gradients
                    dvdξ = 0.0; dvdη = 0.0; dvdζ = 0.0
                    @turbo for ii = 1:ngl
                        dvdξ += dψ[ii,k]*uprimitiveieq[3,ii,l,m]
                        dvdη += dψ[ii,l]*uprimitiveieq[3,k,ii,m]
                        dvdζ += dψ[ii,m]*uprimitiveieq[3,k,l,ii]
                    end

                    # w-velocity gradients (NEW)
                    dwdξ = 0.0; dwdη = 0.0; dwdζ = 0.0
                    @turbo for ii = 1:ngl
                        dwdξ += dψ[ii,k]*uprimitiveieq[4,ii,l,m]
                        dwdη += dψ[ii,l]*uprimitiveieq[4,k,ii,m]
                        dwdζ += dψ[ii,m]*uprimitiveieq[4,k,l,ii]
                    end

                    # Metric terms
                    dξdx_klm = dξdx[k, l, m, iel]
                    dξdy_klm = dξdy[k, l, m, iel]
                    dξdz_klm = dξdz[k, l, m, iel]

                    dηdx_klm = dηdx[k, l, m, iel]
                    dηdy_klm = dηdy[k, l, m, iel]
                    dηdz_klm = dηdz[k, l, m, iel]

                    dζdx_klm = dζdx[k, l, m, iel]
                    dζdy_klm = dζdy[k, l, m, iel]
                    dζdz_klm = dζdz[k, l, m, iel]

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
                                                            uprimitiveieq[1,k,l,m],
                                                            dudx, dvdy, dwdz,
                                                            dudy, dvdx,
                                                            dudz, dwdx,
                                                            dvdz, dwdy,
                                                            0.0,
                                                            0.0,
                                                            PHYS_CONST, Δ2,
                                                            inputs,
                                                            VT, SD)

                        # Stress tensor for u-momentum
                        τ_xx = 2.0 * effective_viscosity * dudx - (2.0/3.0) * effective_viscosity * div_u
                        τ_xy = effective_viscosity * (dudy + dvdx)
                        τ_xz = effective_viscosity * (dudz + dwdx)

                        flux_x = τ_xx
                        flux_y = τ_xy
                        flux_z = τ_xz
                        μ_local = effective_viscosity

                    elseif is_v_momentum
                        # USE EFFECTIVE VISCOSITY
                        effective_viscosity = SGS_diffusion(visc_coeffieq, ieq,
                                                            uprimitiveieq[1,k,l,m],
                                                            dudx, dvdy, dwdz,
                                                            dudy, dvdx,
                                                            dudz, dwdx,
                                                            dvdz, dwdy,
                                                            0.0,
                                                            0.0,
                                                            PHYS_CONST, Δ2,
                                                            inputs,
                                                            VT, SD)

                        # Stress tensor for v-momentum
                        τ_xy = effective_viscosity * (dudy + dvdx)
                        τ_yy = 2.0 * effective_viscosity * dvdy - (2.0/3.0) * effective_viscosity * div_u
                        τ_yz = effective_viscosity * (dvdz + dwdy)

                        flux_x = τ_xy
                        flux_y = τ_yy
                        flux_z = τ_yz
                        μ_local = effective_viscosity

                    elseif is_w_momentum  # NEW BLOCK
                        # USE EFFECTIVE VISCOSITY
                        effective_viscosity = SGS_diffusion(visc_coeffieq, ieq,
                                                            uprimitiveieq[1,k,l,m],
                                                            dudx, dvdy, dwdz,
                                                            dudy, dvdx,
                                                            dudz, dwdx,
                                                            dvdz, dwdy,
                                                            0.0,
                                                            0.0,
                                                            PHYS_CONST, Δ2,
                                                            inputs,
                                                            VT, SD)

                        # Stress tensor for w-momentum
                        τ_xz = effective_viscosity * (dudz + dwdx)
                        τ_yz = effective_viscosity * (dvdz + dwdy)
                        τ_zz = 2.0 * effective_viscosity * dwdz - (2.0/3.0) * effective_viscosity * div_u

                        flux_x = τ_xz
                        flux_y = τ_yz
                        flux_z = τ_zz
                        μ_local = effective_viscosity

                    elseif is_temperature

                        if inputs[:energy_equation] == "theta"
                            # Compute temperature gradient
                            dθdξ = 0.0; dθdη = 0.0; dθdζ = 0.0
                            @turbo for ii = 1:ngl
                                dθdξ += dψ[ii,k]*uprimitiveieq[ieq,ii,l,m]
                                dθdη += dψ[ii,l]*uprimitiveieq[ieq,k,ii,m]
                                dθdζ += dψ[ii,m]*uprimitiveieq[ieq,k,l,ii]
                            end

                            # Transform to physical coordinates
                            dθdx = dθdξ*dξdx_klm + dθdη*dηdx_klm + dθdζ*dζdx_klm
                            dθdy = dθdξ*dξdy_klm + dθdη*dηdy_klm + dθdζ*dζdy_klm
                            dθdz = dθdξ*dξdz_klm + dθdη*dηdz_klm + dθdζ*dζdz_klm

                            if inputs[:lrichardson]
                                θ_ref = uprimitiveieq[5,k,l,m]  # Local temperature
                            else
                                θ_ref = 1.0  # Dummy value (not used when lrichardson=false)
                            end

                            # USE EFFECTIVE DIFFUSIVITY
                            effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                                uprimitiveieq[1,k,l,m],
                                                                dudx, dvdy, dwdz,
                                                                dudy, dvdx,
                                                                dudz, dwdx,
                                                                dvdz, dwdy,
                                                                θ_ref,
                                                                dθdz,
                                                                PHYS_CONST, Δ2,
                                                                inputs,
                                                                VT, SD)
                            flux_x = effective_diffusivity * dθdx
                            flux_y = effective_diffusivity * dθdy
                            flux_z = effective_diffusivity * dθdz
                            μ_local = effective_diffusivity

                        elseif inputs[:energy_equation] == "energy"
                            PhysConst = PhysicalConst{Float32}()
                            cp        = PhysConst.cp
                            Rvap      = PhysConst.Rvap
                            Lc        = PhysConst.Lc
                            ip        = connijk[iel,k,l,m]
                            # Compute energy gradient
                            dhldξ = 0.0; dhldη = 0.0; dhldζ = 0.0
                            @turbo for ii = 1:ngl
                                dhldξ += dψ[ii,k]*uprimitiveieq[ieq,ii,l,m]
                                dhldη += dψ[ii,l]*uprimitiveieq[ieq,k,ii,m]
                                dhldζ += dψ[ii,m]*uprimitiveieq[ieq,k,l,ii]
                            end
                            # Transform to physical coordinates
                            dhldx = dhldξ*dξdx_klm + dhldη*dηdx_klm + dhldζ*dζdx_klm
                            dhldy = dhldξ*dξdy_klm + dhldη*dηdy_klm + dhldζ*dζdy_klm
                            dhldz = dhldξ*dξdz_klm + dhldη*dηdz_klm + dhldζ*dζdz_klm
                            if inputs[:lrichardson]
                                T_ref = Tabs[ip]
                                # θ_ref = Tabs[ip]*(PhysConst.pref/uaux[ip,end])^(1/PhysConst.cpoverR)

                                # Compute condensate mixing ratio gradient
                                dqndξ = 0.0; dqndη = 0.0; dqndζ = 0.0
                                # dθ_refdξ = 0.0; dθ_refdη = 0.0; dθ_refdζ = 0.0
                                # p = uaux[:,end]
                                @turbo for ii = 1:ngl
                                    ip_k  = conn_el[ii,l,m]
                                    ip_l  = conn_el[k,ii,m]
                                    ip_m  = conn_el[k,l,ii]
                                    dqndξ += dψ[ii,k]*qn[ip_k]
                                    dqndη += dψ[ii,l]*qn[ip_l]
                                    dqndζ += dψ[ii,m]*qn[ip_m]
                                    # dθ_refdξ += dψ[ii,k]*Tabs[ip_k]*(PhysConst.pref/p[ip_k])^(1/PhysConst.cpoverR)
                                    # dθ_refdη += dψ[ii,l]*Tabs[ip_l]*(PhysConst.pref/p[ip_l])^(1/PhysConst.cpoverR)
                                    # dθ_refdζ += dψ[ii,m]*Tabs[ip_m]*(PhysConst.pref/p[ip_m])^(1/PhysConst.cpoverR)
                                end
                                # Transform to physical coordinates
                                dqndz = dqndξ*dξdz_klm + dqndη*dηdz_klm + dqndζ*dζdz_klm
                                # dθ_refdz = dθ_refdξ*dξdz_klm + dθ_refdη*dηdz_klm + dθ_refdζ*dζdz_klm

                                γ          = (Lc^2 * qs[ip]) / (Rvap * cp * T_ref^2)
                                dhl_eff_dz =(1.0 / (cp * (1 + γ))) * dhldz - T_ref * dqndz
                            else
                                T_ref      = 1.0 # Dummy value (not used when lrichardson=false)
                                dhl_eff_dz = 1.0
                            end

                             # USE EFFECTIVE DIFFUSIVITY
                            effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                                uprimitiveieq[1,k,l,m],
                                                                dudx, dvdy, dwdz,
                                                                dudy, dvdx,
                                                                dudz, dwdx,
                                                                dvdz, dwdy,
                                                                T_ref,
                                                                dhl_eff_dz,
                                                                PHYS_CONST, Δ2,
                                                                inputs,
                                                                VT, SD)
                            flux_x = effective_diffusivity * dhldx
                            flux_y = effective_diffusivity * dhldy
                            flux_z = effective_diffusivity * dhldz
                            μ_local = effective_diffusivity
                        end


                    else
                        # Other scalars (use appropriate Schmidt number)
                        # USE EFFECTIVE DIFFUSIVITY
                        effective_diffusivity = SGS_diffusion(visc_coeffieq, ieq,
                                                              uprimitiveieq[1,k,l,m],
                                                              dudx, dvdy, dwdz,
                                                              dudy, dvdx,
                                                              dudz, dwdx,
                                                              dvdz, dwdy,
                                                              0.0,
                                                              0.0,
                                                              PHYS_CONST, Δ2,
                                                              inputs,
                                                              VT, SD)

                        # Compute scalar gradient
                        dqdξ = 0.0; dqdη = 0.0; dqdζ = 0.0
                        @turbo for ii = 1:ngl
                            dqdξ += dψ[ii,k]*uprimitiveieq[ieq,ii,l,m]
                            dqdη += dψ[ii,l]*uprimitiveieq[ieq,k,ii,m]
                            dqdζ += dψ[ii,m]*uprimitiveieq[ieq,k,l,ii]
                        end

                        # Transform to physical coordinates
                        dqdx = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
                        dqdy = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
                        dqdz = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm

                        flux_x = effective_diffusivity * dqdx
                        flux_y = effective_diffusivity * dqdy
                        flux_z = effective_diffusivity * dqdz
                        μ_local = effective_diffusivity
                    end

                    # ===== Weak form assembly (3D) =====
                    ∇ξ_flux_klm = (dξdx_klm*flux_x + dξdy_klm*flux_y + dξdz_klm*flux_z)*ωJac * σμ
                    ∇η_flux_klm = (dηdx_klm*flux_x + dηdy_klm*flux_y + dηdz_klm*flux_z)*ωJac * σμ
                    ∇ζ_flux_klm = (dζdx_klm*flux_x + dζdy_klm*flux_y + dζdz_klm*flux_z)*ωJac * σμ

                    @turbo for i = 1:ngl
                        dhdξ_ik = dψ[i,k]
                        dhdη_il = dψ[i,l]
                        dhdζ_im = dψ[i,m]

                        rhs_diffξ_el[iel,i,l,m,ieq] -= dhdξ_ik * ∇ξ_flux_klm
                        rhs_diffη_el[iel,k,i,m,ieq] -= dhdη_il * ∇η_flux_klm
                        rhs_diffζ_el[iel,k,l,i,ieq] -= dhdζ_im * ∇ζ_flux_klm
                    end
                    μ_max_ieq = max(μ_local * σμ, μ_max_ieq)
                end
            end
        end
    end
    μ_max[ieq] = μ_max_ieq
end

function  _expansion_visc!(rhs_diffξ_el, rhs_diffη_el, uprimitiveieq, visc_coeff, ω, mesh, basis, metrics, inputs, rhs_el, iel, ieq, QT::Exact, VT, SD::NSD_2D, ::FD)
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
                @inbounds ωJac = ω[i]*ω[j]*ω[k]*Je[i, j, k, iel]

                dHdξ = 0.0
                dHdη = 0.0
                dHdζ = 0.0
                @turbo for m = 1:ngl
                    dHdξ += dψ[m,i]*q[m,j,k,1]
                    dHdη += dψ[m,j]*q[i,m,k,1]
                    dHdζ += dψ[m,k]*q[i,j,m,1]
                end
                dξdz_ij = dξdz[i, j, k, iel]
                dηdz_ij = dηdz[i, j, k, iel]
                dζdz_ij = dζdz[i, j, k, iel]

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
            ωJac = ω[i]*ω[j]*Je[i, j, iel]

            dHdξ = 0.0
            dHdη = 0.0
            @turbo for m = 1:ngl
                dHdξ += dψ[m,i]*q[m,j,1]
                dHdη += dψ[m,j]*q[i,m,1]
            end
            dξdy_ij = dξdy[i, j, iel]
            dηdy_ij = dηdy[i, j, iel]

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
