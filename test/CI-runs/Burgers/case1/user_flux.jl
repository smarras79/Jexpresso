function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)

    #
    # 1D viscous Burgers equation in conservation form:
    #
    #   ∂q/∂t + ∂F(q)/∂x = ν ∂²q/∂x²,    F(q) = q²/2
    #
    F[1] = 0.5*q[1]*q[1]

end

function user_flux_gpu(q, qe, PhysConst, lpert)

    T = eltype(q)
    return T(0.5*q[1]*q[1]), T(0.0)

end
