function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
    uprimitive[5] = u[5]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = u[4]/(u[1]+qe[1])
    uprimitive[5] = (u[5]+qe[5])/(u[1]+qe[1])-qe[5]/qe[1]
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)

    PhysConst = PhysicalConst{Float64}()
    
    uout[1] = u[1]      #ρ
    uout[2] = u[2]/u[1] #u
    uout[3] = u[3]/u[1] #v
    uout[4] = u[4]/u[1] #w
    uout[5] = u[5]/u[1] #θ
    uout[end] = perfectGasLaw_ρθtoP(PhysConst; ρ=uout[1], θ=uout[5]) #P

end
function user_les_profiles!(means, prof, q, qe, ET)

    PhysConst = PhysicalConst{Float64}()

    if ET == PERT()
        ρ = q[1] + qe[1]
        θ = (q[5] + qe[5]) / ρ - qe[5] / qe[1]
    else
        ρ = q[1]
        θ = q[5] / ρ
    end

    u = q[2] / ρ
    v = q[3] / ρ 
    w = q[4] / ρ 
    p = perfectGasLaw_ρθtoP(PhysConst; ρ=ρ, θ=θ)

    means[1] = u                                          # u
    means[2] = v                                          # v
    means[3] = w                                          # w
    means[4] = θ                                          # θ
    means[5] = p                                          # p

    # ---- resolved Reynolds stress tensor ----
    prof[1]  = u * u  # <uu> res
    prof[2]  = u * v  # <uv> res
    prof[3]  = u * w  # <uw> res
    prof[4]  = v * v  # <vv> res
    prof[5]  = v * w  # <vw> res
    prof[6]  = w * w  # <ww> res

    # ---- resolved heat fluxes and temperature variance ----
    prof[7]  = θ * θ  # <θθ> res
    prof[8]  = u * θ  # <uθ> res
    prof[9]  = v * θ  # <vθ> res
    prof[10] = w * θ  # <wθ> res

    # ---- SFS stress — requires SGS model output; zero until implemented ----
    prof[11] = 0.0  # <uu> sfs
    prof[12] = 0.0  # <uv> sfs
    prof[13] = 0.0  # <uw> sfs
    prof[14] = 0.0  # <vv> sfs
    prof[15] = 0.0  # <vw> sfs
    prof[16] = 0.0  # <ww> sfs
    prof[17] = 0.0  # <θθ> sfs (no SGS scalar variance in Smagorinsky)
    prof[18] = 0.0  # <uθ> sfs
    prof[19] = 0.0  # <vθ> sfs
    prof[20] = 0.0  # <wθ> sfs

    # ---- pressure-velocity correlations ----
    prof[21] = u * p  # <up>
    prof[22] = v * p  # <vp>
    prof[23] = w * p  # <wp>

    # ---- dissipation rates — require velocity/temperature gradients; zero until implemented ----
    prof[24] = 0.0  # eps   (TKE dissipation: 2 ν_eff S_ij S_ij)
    prof[25] = 0.0  # eps_t (θ-variance dissipation: κ_eff |∇θ|²)

    # ---- density ----
    prof[26] = ρ    # rho

    # ---- triple moments ----
    u2 = u * u
    v2 = v * v
    w2 = w * w

    prof[27] = u2 * u  # <uuu>
    prof[28] = u2 * v  # <uuv>
    prof[29] = u2 * w  # <uuw>
    prof[30] = v2 * u  # <vvu>
    prof[31] = v2 * v  # <vvv>
    prof[32] = v2 * w  # <vvw>
    prof[33] = w2 * u  # <wwu>
    prof[34] = w2 * v  # <wwv>
    prof[35] = w2 * w  # <www>
    prof[36] = u2 * θ  # <uuθ>
    prof[37] = v2 * θ  # <vvθ>
    prof[38] = w2 * θ  # <wwθ>

end

function user_les_stress!(profp, prof, means)

    ubar = means[1]  # ū
    vbar = means[2]  # v̄
    wbar = means[3]  # w̄
    θbar = means[4]  # θ̄
    pbar = means[5]  # p̄

    # ---- resolved Reynolds stress tensor ----
    profp[1]  = prof[1] - ubar * ubar  # <u'u'> res
    profp[2]  = prof[2] - ubar * vbar  # <u'v'> res
    profp[3]  = prof[3] - ubar * wbar  # <u'w'> res
    profp[4]  = prof[4] - vbar * vbar  # <v'v'> res
    profp[5]  = prof[5] - vbar * wbar  # <v'w'> res
    profp[6]  = prof[6] - wbar * wbar  # <w'w'> res

    # ---- resolved heat fluxes and temperature variance ----
    profp[7]  = prof[7]  - θbar * θbar  # <θ'θ'> res
    profp[8]  = prof[8]  - ubar * θbar  # <u'θ'> res
    profp[9]  = prof[9]  - vbar * θbar  # <v'θ'> res
    profp[10] = prof[10] - wbar * θbar  # <w'θ'> res

    # ---- SFS stress — requires SGS model output; zero until implemented ----
    profp[11] = prof[11] - 0.0  # <u'u'> sfs
    profp[12] = prof[12] - 0.0  # <u'v'> sfs
    profp[13] = prof[13] - 0.0  # <u'w'> sfs
    profp[14] = prof[14] - 0.0  # <v'v'> sfs
    profp[15] = prof[15] - 0.0  # <v'w'> sfs
    profp[16] = prof[16] - 0.0  # <w'w'> sfs
    profp[17] = prof[17] - 0.0  # <θ'θ'> sfs (no SGS scalar variance in Smagorinsky)
    profp[18] = prof[18] - 0.0  # <u'θ'> sfs
    profp[19] = prof[19] - 0.0  # <v'θ'> sfs
    profp[20] = prof[20] - 0.0  # <w'θ'> sfs

    # ---- pressure-velocity correlations ----
    profp[21] = prof[21] - ubar * pbar  # <u'p'>
    profp[22] = prof[22] - vbar * pbar  # <v'p'>
    profp[23] = prof[23] - wbar * pbar  # <w'p'>

    # ---- dissipation rates — require velocity/temperature gradients; zero until implemented ----
    profp[24] = prof[24] - 0.0  # eps   (TKE dissipation: 2 ν_eff S_ij S_ij)
    profp[25] = prof[25] - 0.0  # eps_t (θ-variance dissipation: κ_eff |∇θ|²)

    # ---- density ----
    profp[26] = prof[26]    # <ρ> — mean density from first pass

    # ---- triple moments ----
    # Identity: <a'a'b'> = <a²b> - b̄<a²> - 2ā<ab> + 2ā²b̄
    # Identity: <a'a'a'> = <a³>  - 3ā<a²> + 2ā³
    ubar2 = ubar * ubar
    vbar2 = vbar * vbar
    wbar2 = wbar * wbar

    profp[27] = prof[27] - 3*ubar*prof[1]  + 2*ubar2*ubar                    # <u'u'u'>
    profp[28] = prof[28] - vbar*prof[1]    - 2*ubar*prof[2] + 2*ubar2*vbar   # <u'u'v'>
    profp[29] = prof[29] - wbar*prof[1]    - 2*ubar*prof[3] + 2*ubar2*wbar   # <u'u'w'>
    profp[30] = prof[30] - ubar*prof[4]    - 2*vbar*prof[2] + 2*vbar2*ubar   # <v'v'u'>
    profp[31] = prof[31] - 3*vbar*prof[4]  + 2*vbar2*vbar                    # <v'v'v'>
    profp[32] = prof[32] - wbar*prof[4]    - 2*vbar*prof[5] + 2*vbar2*wbar   # <v'v'w'>
    profp[33] = prof[33] - ubar*prof[6]    - 2*wbar*prof[3] + 2*wbar2*ubar   # <w'w'u'>
    profp[34] = prof[34] - vbar*prof[6]    - 2*wbar*prof[5] + 2*wbar2*vbar   # <w'w'v'>
    profp[35] = prof[35] - 3*wbar*prof[6]  + 2*wbar2*wbar                    # <w'w'w'>
    profp[36] = prof[36] - θbar*prof[1]    - 2*ubar*prof[8]  + 2*ubar2*θbar  # <u'u'θ'>
    profp[37] = prof[37] - θbar*prof[4]    - 2*vbar*prof[9]  + 2*vbar2*θbar  # <v'v'θ'>
    profp[38] = prof[38] - θbar*prof[6]    - 2*wbar*prof[10] + 2*wbar2*θbar  # <w'w'θ'>

end

function user_les_spectral!(spectra, kappa, u_unif, Ly)
    N_unif = size(u_unif, 1)
    nk     = N_unif ÷ 2 + 1
    N2     = Float64(N_unif * N_unif)

    for ivar in 1:4
        signal = copy(u_unif[:, ivar])
        signal .-= sum(signal) / N_unif   # remove y-mean
        ŝ = fft(signal)
        for ik in 1:nk
            ivar == 1 && (kappa[ik] = (ik-1) / Ly)   # fill kappa once
            fac = (ik == 1 || ik == nk) ? 1.0 : 2.0
            spectra[ik, ivar] = fac * abs2(ŝ[ik]) / N2 * Ly
        end
    end
end