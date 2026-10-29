function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)

    #
    # 1D inviscid Burgers equation in conservation form:
    #
    #   ∂q/∂t + ∂F(q)/∂x = 0,    F(q) = q²/2
    #
    # The only viscosity is the DSGS one (user_inputs.jl), which the viscous
    # operator adds as ∂ₓ(ν_DSGS ∂ₓq).
    #
    F[1] = 0.5*q[1]*q[1]

end

function user_flux_gpu(q, qe, PhysConst, lpert)

    T = eltype(q)
    return T(0.5*q[1]*q[1]), T(0.0)

end
