#
# 1D scalar linear advection:  ∂q/∂t + ∂F(q)/∂x = 0,  F(q) = U q
#
# The advection speed is defined once, in a function, and consumed by BOTH the
# analytic flux and the max-wave-speed hook. Duplicating the literal in two
# places is a silent-failure mode: a mismatch leaves the scheme stable and
# convergent but with the wrong amount of interface dissipation, which shows up
# only as a wrong error constant -- exactly the kind of defect an order study
# passes straight through.
#
advection_speed() = 2.0


function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)

    F[1] = advection_speed() * q[1]

end


function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    return T(advection_speed() * q[1]), T(0.0)
end


#
# Max wave-speed hook for the DG numerical flux (Q3).
#
# Called by surface_rhs_el! (kernel/operators/dg_fluxes.jl) as
#     user_max_wave_speed(qL, @view(qe[ipL,:]), SD, SVT; neqs=neqs)
# so the positional signature and the ::TOTAL dispatch must match exactly.
#
# For scalar linear advection the single characteristic speed is U, so
# λ_max = |U| and the Rusanov flux reduces to exact upwinding for either sign.
#
function user_max_wave_speed(q, qe, SD::NSD_1D, ::TOTAL; neqs=1)
    return abs(advection_speed())
end
