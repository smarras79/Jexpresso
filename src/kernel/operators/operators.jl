function interpolate!(qh, qj, ψ)
    #
    # in/out: qh(x)  of size [1:Nx] where Nx is the number of interpolation points
    # in:     qj(ξ)  of size [1:Nξ] where Nξ = polynomial order + 1 (e.g. LGL points)
    # in:     ψ(ξ,x) Lagrange basis of size [1:Nξ, 1:Nx]
    #
    # Interpolation rule using Lagrange basis ψ with zeros at the ξ nodes:
    # q(x)     = ∑ⱼ{1,Nξ} ψⱼ(x)qⱼ
    # ∂q(x)/∂ξ = ∑ⱼ{1,Nξ} ∂ψⱼ(x)/∂x qⱼ
    #
    for ix = 1:length(qh)
        qh[ix] = 0
        for jlgl = 1:length(qj)
            qh[ix] = qh[ix] + ψ[jlgl, ix]*qj[jlgl]
        end
    end
    
    return qh
end
