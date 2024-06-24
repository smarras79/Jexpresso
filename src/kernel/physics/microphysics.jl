function kessler_nocolumn!(mp::St_Microphysics, params,
                           physConst,
                           ρ, t, qv, qc, qr, p, q, z,
                           rainnc, rainncv, dt, nz, qref, mphconst)
    
    # Constants
    c1     = 0.001
    c2     = 0.001
    c3     = 2.2
    c4     = 0.875
    fudge  = 1.0
    mxfall = 10.0

    # Constants related to saturation water pressure
    #mphconst = MicrophysicalConst{Float64}()
    xlv    = mphconst.xlv
    ep2    = mphconst.ep2
    svp1   = mphconst.svp1
    svp2   = mphconst.svp2
    svp3   = mphconst.svp3
    svpt0  = mphconst.svpt0
    ρwater = mphconst.ρwater

    # Global variables
    nvar  = params.neqs
    npoin = params.mesh.npoin
    ngl   = params.mesh.ngl
    nop   = params.mesh.nop

    # Local variables
    qrprod, ern, gam, rcgs, rcgsi = 0.0, 0.0, 0.0, 0.0, 0.0
    prod   = mp.prod
    rhs    = mp.rhs
    vt     = mp.vt
    prodk  = mp.prodk
    vtden  = mp.vtden
    rdzk   = mp.rdzk
    ρk     = mp.ρk
    rdzw   = mp.rdzw
    
    nfall, n, nfall_new, icol                    = 0, 0, 0, 0
    qrr, pressure, temp, es, qvs, dz, pii        = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    f5, dtfall, rdz, product                     = 0.0, 0.0, 0.0, 0.0
    max_heating, max_condense, max_rain, maxqrp  = 0.0, 0.0, 0.0, -100.0
    vtmax, ernmax, crmax, factorn, time_sediment = 0.0, 0.0, 0.0, 0.0, dt
    qcr, factorr, ppt = 0.0, 0.0, 0.0

    # Terminal velocity calculation and advection
    f5 = svp2 * (svpt0 - svp3) * xlv / cp
    max_heating, max_condense, max_rain = 0.0, 0.0, 0.0
    crmax = 0.0

    # Sedimentation
    for ip = 1:npoin
        rdzk[ip]  = 1.0 / ((ztop - zbottom) / max(nelz * nop, 1))
        prodk[ip] = qr[ip]
        ρk[ip]    = ρ[ip]
        qrr       = max(0.0, qr[ip] * 0.001 * ρk[ip])
        vtden[ip] = sqrt(1.123 / ρk[ip])
        vt[ip]    = 36.34 * (qrr^0.1364) * vtden[ip]
        crmax     = max(vt[ip] * dt * rdzk[ip], crmax)
    end

    nfall = max(1, round(Int, 0.5 + crmax / max_cr_sedimentation))
    dtfall = dt / nfall
    time_sediment = dt

    # Terminal velocity calculation and advection loop
    q0 = copy(prodk)
    column_sedimentation = true
    while nfall > 0
        time_sediment -= dtfall
        for ip = 1:npoin
            ppt = 0.0
            ppt = ρk[ip] * prodk[ip] * vt[ip] * dtfall / ρwater
            if coord[3, ip] < 1.0
                rainncv[ip] = ppt * 1000.0
                rainnc[ip] += ppt * 1000.0
            end
        end
        # Call to ti_rk35_production_no_column not translated as its definition is not provided
        if nfall > 1
            nfall -= 1
            crmax = 0.0
            for ip = 1:npoin
                qrr = max(0.0, prodk[ip] * 0.001 * ρk[ip])
                vt[ip] = 36.34 * (qrr^0.1346) * vtden[ip]
                crmax = max(vt[ip] * time_sediment * rdzw[ip], crmax)
            end
            nfall_new = max(1, round(Int, 0.5 + crmax / max_cr_sedimentation))
            if nfall_new != nfall
                nfall = nfall_new
                dtfall = time_sediment / nfall
            end
        else
            for ip = 1:npoin
                prod[ip] = prodk[ip]
            end
            nfall = 0
        end
        q0 = copy(prodk)
    end

    # Conversion processes
    for ip = 1:npoin
        factorn = 1.0 / (1.0 + c3 * dt * max(0.0, qr[ip])^c4)
        qrprod = qc[ip] * (1.0 - factorn) + factorn * c1 * dt * max(qc[ip] - c2, 0.0)
        rcgs = 0.001 * ρ[ip]
        qc[ip] = max(qc[ip] - qrprod, 0.0)
        qr[ip] = qr[ip] + prod[ip] - qr[ip]
        qr[ip] = max(qr[ip] + qrprod, 0.0)
        pii = (p[ip] / p00)^(rgas / cp)
        temp = t[ip] * pii
        pressure = p[ip]
        gam = 2.5e6 / (cp * pii)
        es = 1000.0 * svp1 * exp(svp2 * (temp - svpt0) / (temp - svp3))
        qvs = ep2 * es / (pressure - es)
        prod[ip] = (qv[ip] - qvs) / (1.0 + pressure / (pressure - es) * qvs * f5 / (temp - svp3)^2)
        ern = min(dt * ((1.6 + 124.9 * (rcgs * qr[ip])^0.2046) * (rcgs * qr[ip])^0.525) / (2.55e8 / (pressure * qvs) + 5.4e5) * (dim(qvs, qv[ip]) / (rcgs * qvs)), max(-prod[ip] - qc[ip], 0.0), qr[ip])
        product = max(prod[ip], -qc[ip])
        t[ip] += gam * (product - ern)
        qv[ip] = max(qv[ip] - product + ern, 0.0)
        qc[ip] += product
        qr[ip] -= ern
    end
end
