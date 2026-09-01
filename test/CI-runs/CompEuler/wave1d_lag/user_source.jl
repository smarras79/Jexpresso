function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1, x=0.0, y=0.0, xmin=0.0, xmax=1.0)
    
    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #

    if (x >= 2.49)#nsponge_points * dsy) #&& dbl >= 0.0)
        sponge_coe = 2.0/(1+exp((0.3*(xmax-2.5)-x+2.5)/((xmax)/18)))
    elseif (x <=-2.49)
        sponge_coe = 2.0/(1+exp(-(0.3*(xmin+2.5)-x-2.5)/((xmax)/18)))
    else
        sponge_coe = 0.0
    end
    #
    # SPONGE STRENGTH: the 20 is not a retune, it is the bug fix bookkeeping.
    #
    # The 1D SEM kernels used to leave the element Jacobian off the source
    # term while the mass-matrix division still carried it, so every 1D source
    # actually reached the ODE as S/Je. Here Je = 0.05 BOTH in the CG region
    # (Δx/2 = (2.5 - (-2.5))/50/2) and in the Laguerre region (Je is set to
    # :yfac_laguerre = 0.05), so this sponge ran at exactly 1/0.05 = 20× the
    # coefficient written below, everywhere it is active.
    #
    # The kernels now apply Je correctly — see the note above
    # _expansion_inviscid!(..., NSD_1D, ContGal) in
    # src/kernel/operators/rhs.jl — so the 20 is folded in here to hold this
    # case's absorbing layer at the strength it was tuned to. Because a single
    # Je covers both regions, the compensation is exact and the solution is
    # unchanged to round-off.
    #
    cs = min(sponge_coe,1)/0.05
    S[1] = -(cs)*(q[1])
    S[2] = -(cs)*(q[2])
    return  S
end

function user_source_gpu(q,qe,x,PhysConst, xmax, xmin,lpert)

    T = eltype(q)
    if (x >= 2.49)#nsponge_points * dsy) #&& dbl >= 0.0)
        sponge_coe = 2.0/(1+exp((0.3*(xmax-2.5)-x+2.5)/((xmax)/18)))
    elseif (x <=-2.49)
        sponge_coe = 2.0/(1+exp(-(0.3*(xmin+2.5)-x-2.5)/((xmax)/18)))
    else
        sponge_coe = 0.0
    end
    cs = min(sponge_coe,1)/T(0.05)   # see the note in user_source! above
    return T(-cs*q[1]), T(-cs*q[2])

end
