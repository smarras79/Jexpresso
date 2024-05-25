function kessler(mp::St_Microphysics, params, q, qref, mesh)
    
    temp3 = 0.0
    temp4 = 0.0
    
    #= allocate if not already allocated
    if isempty(rho)
        rho = zeros(Float64, npoin)
        t = zeros(Float64, npoin)
        p = zeros(Float64, npoin)
        z = zeros(Float64, npoin)
        rnc = zeros(Float64, npoin)
        rncv = zeros(Float64, npoin)
        temp1 = zeros(Float64, npoin)
        temp2 = zeros(Float64, npoin)
    end=#
    
    #temp1 .= rnc
    #temp2 .= rncv
    
    kessler_nocolumn!(mp,
                      params,
                      mesh.npoin,
                      @view(q[:,1]),     #ρ
                      @view(q[:,mesh.nsd+2]), #θ
                      @view(q[:,mesh.nsd+2]), #qv
                      @view(q[:,mesh.nsd+3]), #qc
                      @view(q[:,mesh.nsd+4]), #qr
                      qref,
                      mesh.z,
                      params.Δt)

    #for ip in 1:npoin
    #    # Store accumulated rain
    #    rainnc[ip] = temp1[ip]
    #    rainncv[ip] = temp2[ip]
    #end
    
end

function kessler_nocolumn!(mp::St_Microphysics,
                           params,
                           npoin,
                           ρ,
                           t,
                           qv,
                           qc,
                           qr,
                           qref,
                           z,
                           dt)
    
    # Constants
    c1     = 0.001
    c2     = 0.001
    c3     = 2.2
    c4     = 0.875
    fudge  = 1.0
    mxfall = 10.0

    # Constants related to saturation water pressure
    physConst = PhysicalConst{Float64}()
    cp        = physConst.cp
    
    mphConst  = MicrophysicalConst{Float64}()
    xlv       = mphConst.xlv
    ep2       = mphConst.ep2
    svp1      = mphConst.svp1
    svp2      = mphConst.svp2
    svp3      = mphConst.svp3
    svpt0     = mphConst.svpt0
    ρwater    = mphConst.ρwater

    # Global variables
    nvar    = params.neqs
    npoin   = params.mesh.npoin
    ngl     = params.mesh.ngl
    nop     = params.mesh.nop
    zbottom = params.zmin
    ztop    = params.zmax

    # Local variables
    qrprod, ern, gam, rcgs, rcgsi = 0.0, 0.0, 0.0, 0.0, 0.0
    prod   = mp.prod
    prodk  = mp.prodk
    vt     = mp.vt
    vtden  = mp.vtden
    rdzk   = mp.rdzk
    rdzw   = mp.rdzw
    
    nfall, n, nfall_new, icol                    = 0, 0, 0, 0
    qrr, pressure, temp, es, qvs, dz, pii        = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    f5, dtfall, rdz, product                     = 0.0, 0.0, 0.0, 0.0
    max_heating, max_condense, max_rain, maxqrp  = 0.0, 0.0, 0.0, -100.0
    vtmax, ernmax, crmax, factorn                = 0.0, 0.0, 0.0, 0.0
    time_sediment = dt
    qcr, factorr, ppt = 0.0, 0.0, 0.0
    max_cr_sedimentation = 0.75
    
    # Terminal velocity calculation and advection
    f5           = svp2 * (svpt0 - svp3) * xlv / cp
    max_heating  = 0.0
    max_condense = 0.0
    max_rain     = 0.0
    crmax        = 0.0

    # Sedimentation
    for ip = 1:npoin
        rdzk[ip]  = 1.0 / ((ztop - zbottom) / max(nelz * nop, 1))
        prodk[ip] = qr[ip]
        qrr       = max(0.0, qr[ip] * 0.001 * ρ[ip])
        vtden[ip] = sqrt(1.123 / ρ[ip])
        vt[ip]    = 36.34 * (qrr^0.1364) * vtden[ip]
        crmax     = max(vt[ip] * dt * rdzk[ip], crmax)
    end
    
#=
    nfall         = max(1, Int64(floor((0.5 + crmax / max_cr_sedimentation))))
    dtfall        = dt / nfall
    time_sediment = dt

    # Terminal velocity calculation and advection loop
    q0 = copy(prodk)
    column_sedimentation = true
    while nfall > 0
        time_sediment -= dtfall
        for ip = 1:npoin
            ppt = 0.0
            ppt = ρ[ip] * prodk[ip] * vt[ip] * dtfall / ρwater
            #if coord[3, ip] < 1.0
            #    rainncv[ip] = ppt * 1000.0
            #    rainnc[ip] += ppt * 1000.0
            #end
        end
        # Call to ti_rk35_production_no_column not translated as its definition is not provided
        if nfall > 1
            nfall -= 1
            crmax = 0.0
            for ip = 1:npoin
                qrr = max(0.0, prodk[ip] * 0.001 * ρ[ip])
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
    =#
end
