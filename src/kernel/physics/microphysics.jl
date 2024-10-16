function do_micro_physics!(mp::St_SamMicrophysics, npoin, uaux, z, qe, ::TOTAL)

     PhysConst = PhysicalConst{TFloat}()
     MicroConst = MicrophysicalConst{TFloat}()

    #get pressure and perform saturation adjustment
    for ip=1:npoin
        uaux[ip,6] = max(0.0,uaux[ip,6])
        uaux[ip,7] = max(0.0,uaux[ip,7])
        T,P,qn,qc,qi,qr,qs,qg,qsatt = saturation_adjustment_sam_microphysics(uaux[ip,1],uaux[ip,5]/uaux[ip,1],uaux[ip,6]/uaux[ip,1],uaux[ip,7]/uaux[ip,1],
                                                                             z[ip],MicroConst,PhysConst)
        mp.Tabs[ip] = T
        mp.qn[ip] = qn
        mp.qi[ip] = qi
        mp.qc[ip] = qc
        mp.qr[ip] = qr
        mp.qs[ip] = qs
        mp.qg[ip] = qg
        uaux[ip,5] = uaux[ip,1]*(PhysConst.cp*T + PhysConst.g*z[ip] - MicroConst.Lc*(qc + qr) - MicroConst.Ls*(qi + qs + qg))
        uaux[ip,end] = P

        ### find Pr, Ps, Pg
        Pr, Ps, Pg = compute_Pm(uaux[ip,1],qr, qs, qg, MicroConst)
        mp.Pr[ip] = Pr
        mp.Ps[ip] = Ps
        mp.Pg[ip] = Pg
        
        ### compute dqp/dt from microphysical processes
        
        mp.S_micro[ip] = compute_dqpdt_sam_micro(uaux[ip,1],T,P,uaux[ip,5]/uaux[ip,1],uaux[ip,6]/uaux[ip,1],qn,qc,qi,uaux[ip,7]/uaux[ip,1],qr,qs,qg,qsatt,MicroConst,PhysConst)

    end

    ### find vertical derivatives for sedimentation processes. 
    
    return nothing
end

function do_micro_physics!(mp::St_SamMicrophysics, npoin, uaux, z, qe, ::PERT)

     PhysConst = PhysicalConst{TFloat}()
     MicroConst = MicrophysicalConst{TFloat}()

    #get pressure and perform saturation adjustment
    for ip=1:npoin
        ρqt_total = max(0.0,uaux[ip,6]+qe[ip,6])
        ρqp_total = max(0.0,uaux[ip,7]+qe[ip,7])
        uaux[ip,6] = ρqt_total - qe[ip,6]
        uaux[ip,7] = ρqp_total - qe[ip,7]
        T,P,qn,qc,qi,qr,qs,qg,qsatt = saturation_adjustment_sam_microphysics(uaux[ip,1]+qe[ip,1],(uaux[ip,5]+qe[ip,5])/(uaux[ip,1]+qe[ip,1]),
                                                                             (uaux[ip,6]+qe[ip,6])/(uaux[ip,1]+qe[ip,1]),
                                                                             (uaux[ip,7]+qe[ip,7])/(uaux[ip,1]+qe[ip,1]),
                                                                             z[ip],MicroConst,PhysConst)
        mp.Tabs[ip] = T
        mp.qn[ip] = qn
        mp.qi[ip] = qi
        mp.qc[ip] = qc
        mp.qr[ip] = qr
        mp.qs[ip] = qs
        mp.qg[ip] = qg
        #hl = (uaux[ip,1]+qe[ip,1])*(PhysConst.cp*T + PhysConst.g*z[ip] - MicroConst.Lc*(qc + qr) - MicroConst.Ls*(qi+qs+qg))
        #uaux[ip,5] = hl - qe[ip,5]
        ## retrieve absolute temp
        #T = (hl/(uaux[ip,1]+qe[ip,1]) - PhysConst.g*z[ip] + MicroConst.Lc*(qc + qr) + MicroConst.Ls*(qi + qs + qg))/PhysConst.cp 
        #mp.Tabs[ip] = T
        #P = moistPressure(PhysConst; ρ = uaux[ip,1]+qe[ip,1], Temp = T, qv = (uaux[ip,6]+qe[ip,6])/(uaux[ip,1]+qe[ip,1])-qn)
        #hl = (uaux[ip,5]+qe[ip,5])/(uaux[ip,1]+qe[ip,1])
        #T = (hl - PhysConst.g*z[ip] + MicroConst.Lc*(qc + qr) + MicroConst.Ls*(qi+qs+qg))/PhysConst.cp
        #P = moistPressure(PhysConst; ρ = uaux[ip,1]+qe[ip,1], Temp = T, qv = (uaux[ip,6]+qe[ip,6])/(uaux[ip,1]+qe[ip,1])-qn)
        #mp.Tabs[ip] = T
        uaux[ip,end] = P

        ### find Pr, Ps, Pg
        Pr, Ps, Pg = compute_Pm(uaux[ip,1]+qe[ip,1], qr, qs, qg, MicroConst)
        mp.Pr[ip] = Pr
        mp.Ps[ip] = Ps
        mp.Pg[ip] = Pg

        ### compute dqp/dt from microphysical processes

        mp.S_micro[ip] = compute_dqpdt_sam_micro(uaux[ip,1]+qe[ip,1],T,P,(uaux[ip,5]+qe[ip,5])/(uaux[ip,1]+qe[ip,1]),
                                                 (uaux[ip,6]+qe[ip,6])/(uaux[ip,1]+qe[ip,1]),qn,qc,qi,
                                                 (uaux[ip,7]+qe[ip,7])/(uaux[ip,1]+qe[ip,1]),qr,qs,qg,qsatt,MicroConst,PhysConst)

    end

    ### find vertical derivatives for sedimentation processes.

    return nothing
end

function compute_precipitation_derivatives(mp::St_SamMicrophysics,ρ,ρe,hl,nelem,ngl,connijk,H,metrics,ω,dψ,::TOTAL)

    MicroConst = MicrophysicalConst{TFloat}()
    mp.dqpdt .= 0.0
    mp.dqtdt .= 0.0
    mp.dhldt .= 0.0
    for e=1:nelem
        # do precipitation
        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    H[i,j,k,1] = mp.Pr[ip] + mp.Ps[ip] + mp.Pg[ip]
                    #if H[i,j,k,1] > 0
                    #    @info mp.Pr[ip], mp.Ps[ip] , mp.Pg[ip]
                    #end
                end
            end
        end
        compute_vertical_derivative_q!(mp.dqpdt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ)
        
        # do precipitation effects on hl

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    H[i,j,k,1] = MicroConst.Lc*mp.Pr[ip] + MicroConst.Ls*(mp.Ps[ip] + mp.Pg[ip])
                end
            end
        end
        compute_vertical_derivative_q!(mp.dhldt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ)
        #do cloud ice sedimentation

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    H[i,j,k,1] = mp.qi[ip] * ρ[ip] * 0.4 
                end
            end
        end
        compute_vertical_derivative_q!(mp.dqtdt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ)
        #@info maximum(mp.dqtdt), maximum(mp.dhldt), maximum(mp.dqpdt)
        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    T = mp.Tabs[ip]
                    ωn = max(0,min(1,(T-MicroConst.T00n)/(MicroConst.T0n - MicroConst.T00n)))
                    mp.dhldt[e,i,j,k] += (MicroConst.Lc + ωn*MicroConst.Lf)*mp.dqtdt[e,i,j,k]
                end
            end
        end
    end

    H[:,:,:,:] .= 0.0
end

function compute_precipitation_derivatives(mp::St_SamMicrophysics,ρ,ρe,hl,nelem,ngl,connijk,H,metrics,ω,dψ,::PERT)

    MicroConst = MicrophysicalConst{TFloat}()
    mp.dqpdt .= 0.0
    mp.dqtdt .= 0.0
    mp.dhldt .= 0.0
    for e=1:nelem
        # do precipitation
        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    H[i,j,k,1] = mp.Pr[ip] + mp.Ps[ip] + mp.Pg[ip]
                    #if H[i,j,k,1] > 0
                    #    @info mp.Pr[ip], mp.Ps[ip] , mp.Pg[ip]
                    #end
                end
            end
        end
        compute_vertical_derivative_q!(mp.dqpdt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ)

        # do precipitation effects on hl

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    H[i,j,k,1] = MicroConst.Lc*mp.Pr[ip] + MicroConst.Ls*(mp.Ps[ip] + mp.Pg[ip])
                end
            end
        end
        compute_vertical_derivative_q!(mp.dhldt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ)
        #do cloud ice sedimentation

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    H[i,j,k,1] = mp.qi[ip] * (ρ[ip] +ρe[ip])* 0.4
                end
            end
        end
        compute_vertical_derivative_q!(mp.dqtdt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ)

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    T = mp.Tabs[ip]
                    ωn = max(0,min(1,(T-MicroConst.T00n)/(MicroConst.T0n - MicroConst.T00n)))
                    mp.dhldt[e,i,j,k] += (MicroConst.Lc + ωn*MicroConst.Lf)*mp.dqtdt[e,i,j,k]
                end
            end
        end
    end

    H[:,:,:,:] .= 0.0
end

function add_micro_precip_sources!(mp::St_SamMicrophysics,T,S_micro,S,q,qn,qe,::TOTAL)

    PhysConst = PhysicalConst{TFloat}()
    ρ = q[1]
    qt = q[6]/ρ
    qp = q[7]/ρ
    qv = qt - qn
    S[4] += -ρ*PhysConst.g*(0.608*qv-qn-qp) #moisture buoyancy contribution
    S[6] += -ρ*S_micro
    S[7] += ρ*S_micro

end

function add_micro_precip_sources!(mp::St_SamMicrophysics,T,S_micro,S,q,qn,qe,::PERT)

    PhysConst = PhysicalConst{TFloat}()
    ρ = q[1]+qe[1]
    qt = (q[6] + qe[6])/ρ 
    qp = q[7]/ρ
    qv = (qt - qn) - qe[6]/qe[1]
    #@info S[4], -q[1]*PhysConst.g*(0.608*qv-qn-qp)
    #@info S[4], S[4]-q[1]*PhysConst.g*(0.608*qv-qn-qp), q[end] - qe[end]
    S[4] += -q[1]*PhysConst.g*(0.608*qv-qn-qp) #moisture buoyancy contribution
    S[6] += -ρ*S_micro
    S[7] += ρ*S_micro

end



function compute_dqpdt_sam_micro(ρ,T,P,hl,qt,qn,qc,qi,qp,qr,qs,qg,qsatt,MicroConst,PhysConst)
    a_rain = MicroConst.a_rain #Constant in fall speed for rain
    b_rain = MicroConst.b_rain
    N0_rain = MicroConst.N0_rain
    Er_c = MicroConst.Er_c
    Er_i = MicroConst.Er_i
    γ3br = MicroConst.γ3br
    ρ_rain = MicroConst.ρ_rain
    C_rain = MicroConst.C_rain
    a_fr   = MicroConst.a_fr
    b_fr   = MicroConst.b_fr
    γ5br   = MicroConst.γ5br

    a_snow = MicroConst.a_snow #Constant in fall speed for rain
    b_snow = MicroConst.b_snow
    N0_snow = MicroConst.N0_snow
    Es_c = MicroConst.Es_c
    Es_i = MicroConst.Es_i
    γ3bs = MicroConst.γ3bs
    ρ_snow = MicroConst.ρ_snow
    C_snow = MicroConst.C_snow
    a_fs   = MicroConst.a_fs
    b_fs   = MicroConst.b_fs
    γ5bs   = MicroConst.γ5bs

    a_graupel = MicroConst.a_graupel #Constant in fall speed for rain
    b_graupel = MicroConst.b_graupel
    N0_graupel = MicroConst.N0_graupel
    Eg_c = MicroConst.Eg_c
    Eg_i = MicroConst.Eg_i
    γ3bg = MicroConst.γ3bg
    ρ_graupel = MicroConst.ρ_graupel
    C_graupel = MicroConst.C_graupel
    a_fg   = MicroConst.a_fg
    b_fg   = MicroConst.b_fg
    γ5bg   = MicroConst.γ5bg

    ρ0 = MicroConst.ρ0
    α = MicroConst.α
    qc0 = MicroConst.qc0
    β = MicroConst.β
    qi0 = MicroConst.qi0

    Lc = MicroConst.Lc
    Ls = MicroConst.Ls
    Ka = MicroConst.Ka
    Da = MicroConst.Da
    Rvap = PhysConst.Rvap
    Rair = PhysConst.Rair
    μ = MicroConst.μ
    e_satw = esatw(T)
    e_sati = esati(T)

    ### Collection of condensates
    Ar_c = π/4*a_rain*N0_rain*Er_c*γ3br*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_rain*N0_rain))^((3+b_rain)/4)
    Ar_i = exp(0.025*(T-273.16))*π/4*a_rain*N0_rain*Er_i*γ3br*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_rain*N0_rain))^((3+b_rain)/4)
    As_c = π/4*a_snow*N0_snow*Es_c*γ3bs*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_snow*N0_snow))^((3+b_snow)/4)
    As_i = exp(0.025*(T-273.16))*π/4*a_snow*N0_snow*Es_i*γ3bs*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_snow*N0_snow))^((3+b_snow)/4)
    Ag_c = π/4*a_graupel*N0_graupel*Eg_c*γ3bs*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_graupel*N0_graupel))^((3+b_graupel)/4)
    Ag_i = exp(0.025*(T-273.16))*π/4*a_graupel*N0_graupel*Eg_c*γ3bs*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_graupel*N0_graupel))^((3+b_graupel)/4)

    dqrdt = (Ar_c*qc + Ar_i*qi)*qr^((3+b_rain)/4)
    dqsdt = (As_c*qc + As_i*qi)*qs^((3+b_snow)/4)
    dqgdt = (Ag_c*qc + Ag_i*qi)*qg^((3+b_graupel)/4)

    ### Autoconversion
    
    Auto = max(0.0, α*(qc - qc0))

    ### Aggregation

    Aggr = max(0.0, β*exp(0.025*(T-273.16))*(qi-qi0))

    ### Evaporation

    S = (qt-qn)/qsatt
    A_rain = (Lc/(Ka*T))*((Lc/(Rvap*T)) - 1) 
    A_snow_graupel = (Ls/(Ka*T))*((Ls/(Rvap*T)) - 1)
    B_rain = Rvap*Rair/(Da*e_satw)
    B_snow_graupel = Rvap*Rair/(Da*e_sati)
    A_er = a_fr*(ρ/(π*ρ_rain*N0_rain))^(0.5)
    A_es = a_fs*(ρ/(π*ρ_snow*N0_snow))^(0.5)
    A_eg = a_fg*(ρ/(π*ρ_graupel*N0_graupel))^(0.5)
    B_er = b_fr*(ρ*a_rain/μ)^(0.5)*γ5br*(ρ0/ρ)^(0.25)*(ρ/(π*ρ_rain*N0_rain))^((5+b_rain)/8)
    B_es = b_fs*(ρ*a_snow/μ)^(0.5)*γ5bs*(ρ0/ρ)^(0.25)*(ρ/(π*ρ_snow*N0_snow))^((5+b_snow)/8)
    B_eg = b_fg*(ρ*a_graupel/μ)^(0.5)*γ5bg*(ρ0/ρ)^(0.25)*(ρ/(π*ρ_graupel*N0_graupel))^((5+b_graupel)/8)
    
    Evap_r = (2*π*C_rain*N0_rain)/(ρ*(A_rain+B_rain))*(A_er*sqrt(qr) + B_er*qr^((5+b_rain)/8))*(S-1)
    Evap_s = (2*π*C_snow*N0_snow)/(ρ*(A_snow_graupel+B_snow_graupel))*(A_es*sqrt(qs) + B_es*qs^((5+b_snow)/8))*(S-1)
    Evap_g = (2*π*C_graupel*N0_graupel)/(ρ*(A_snow_graupel+B_snow_graupel))*(A_eg*sqrt(qg) + B_eg*qg^((5+b_graupel)/8))*(S-1)

    ### cloud ice sedimentation is handled as part of the precipitation package, this routine ignores it
    
    dqpdt = clamp(dqrdt + dqsdt + dqgdt + Auto + Aggr + Evap_r + Evap_s + Evap_g,-qp,qn)

    #if (qn > 0.001)
        #@info dqrdt, dqsdt, dqgdt, Auto, Aggr, Evap_r, Evap_s, Evap_g, qn, qc, qi, qp, qr, qs, qg, dqpdt
    #end

    return dqpdt 
end

function compute_Pm(ρ,qr, qs, qg, MicroConst)
    a_rain = MicroConst.a_rain #Constant in fall speed for rain
    γ4br = MicroConst.γ4br
    ρ_rain = MicroConst.ρ_rain
    N0_rain = MicroConst.N0_rain
    b_rain = MicroConst.b_rain
    ρ0 = MicroConst.ρ0

    a_snow = MicroConst.a_snow #Constant in fall speed for rain
    γ4bs = MicroConst.γ4bs
    ρ_snow = MicroConst.ρ_snow
    N0_snow = MicroConst.N0_snow
    b_snow = MicroConst.b_snow

    a_graupel = MicroConst.a_graupel #Constant in fall speed for rain
    γ4bg = MicroConst.γ4bg
    ρ_graupel = MicroConst.ρ_graupel
    N0_graupel = MicroConst.N0_graupel
    b_graupel = MicroConst.b_graupel

    Pr = (a_rain * γ4br)/6 * (π * ρ_rain * N0_rain)^(-b_rain/4)*(ρ0/ρ)^(0.5)*(ρ*qr)^(1+b_rain/4)    
    Ps = (a_snow * γ4bs)/6 * (π * ρ_snow * N0_snow)^(-b_snow/4)*(ρ0/ρ)^(0.5)*(ρ*qs)^(1+b_snow/4)
    Pg = (a_graupel * γ4bg)/6 * (π * ρ_graupel * N0_graupel)^(-b_graupel/4)*(ρ0/ρ)^(0.5)*(ρ*qg)^(1+b_graupel/4)
    #if (Pr > 0)
    #    @info Pr, Ps, Pg
    #end
    return Pr, Ps, Pg
end

function saturation_adjustment_sam_microphysics(ρ,hl,qt,qp,z,MicroConst,PhysConst)

    T00n = MicroConst.T00n #Temperature threshold for cloud water
    T0n = MicroConst.T0n #Temperature threshold for ice
    T00p = MicroConst.T00p #Temperature threshold for rain     
    T0p = MicroConst.T0p #Temperature threshold for snow/graupel 
    T00g = MicroConst.T00g #Temperature threshold for graupel
    T0g = MicroConst.T0g #Temperature threshold for graupel
    Lc = MicroConst.Lc
    Lf = MicroConst.Lf
    Ls = MicroConst.Ls
    g=PhysConst.g
    cp = PhysConst.cp
    fac_cond = Lc/cp
    fac_fus = Lf/cp
    fac_sub = Ls/cp
    # find equilibrium temperature from saturation adjustment
    # initial guess for sensible temperature and pressure assumes no condensates/all vapor
    T =  (hl - g*z)/cp
    qt = max(0.0,qt)
    qp = max(0.0,qp)
    #P = moistPressure(PhysConst; ρ=ρ, Temp=T, qv = qt) 
    an = 1/(T0n - T00n)
    bn = T00n * an
    ap = 1/(T0p - T00p)
    bp = T00p*ap
    fac1 = fac_cond+(1+bp)*fac_fus
    fac2 = fac_fus*ap
    ag = 1/(T0g - T00g)
    T1 = T + fac1*qp/(1+fac2*qp)
    P = moistPressure(PhysConst; ρ=ρ, Temp=T1, qv = qt)
    #if (qp > 1e-8) 
    #    @info qp, T, T1, fac1*qp/(1+fac2*qp), fac1*qp, (1+fac2*qp)
    #end

    if (T1 >= T0n)

        T1 = T + fac_cond*qp
        P = moistPressure(PhysConst; ρ=ρ, Temp=T1, qv = qt)
        qsatt = qsatw(T1, P/100)

    elseif (T1 <= T00n)

        T1 = T + fac_sub*qp
        P = moistPressure(PhysConst; ρ=ρ, Temp=T1, qv = qt)
        qsatt = qsati(T1, P/100)

    else

        ωn = max(0,min(1,an*T1-bn))
        qsatt = ωn*qsatw(T1,P/100)+(1-ωn)*qsati(T1,P/100)

    end

    if (qt > qsatt)
        
        niter = 0
        dT = 100
        dqsat = 0.0
        while (abs(dT) > 0.001 && niter < 50)
            
            if (T1 >= T0n)
                
                ωn=1
                lstarn = fac_cond
                dlstarn = 0
                qsatt = qsatw(T1,P/100)
                dqsat = dtqsatw(T1,P/100)

            elseif (T1 <= T00n)
                
                ωn = 0
                lstarn = fac_sub
                dlstarn = 0
                qsatt = qsati(T1,P/100)
                dqsat = dtqsati(T1,P/100)

            else

                ωn = max(0,min(1,an*T1-bn))
                lstarn = fac_cond+(1-ωn)*fac_fus
                dlstarn = an*fac_fus
                qsatt = ωn*qsatw(T1,P/100) + (1-ωn)*qsati(T1,P/100)
                dqsat = ωn*dtqsati(T1,P/100) + (1-ωn)*dtqsati(T1,P/100)
            
            end
            
            if (T1 >= T0p)
                
                ωp = 1
                lstarp = fac_cond
                dlstarp = 0
                
            elseif (T1 <= T00p)

                ωp = 0
                lstarp = fac_sub
                dlstarp = 0

            else

                ωp = max(0,min(1,ap*T1-bp))
                lstarp = fac_cond + (1-ωp)*fac_fus
                dlstarp=ap*fac_fus

            end

            fff = T - T1 + lstarn*(qt - qsatt) + lstarp*qp
            dfff = dlstarn*(qt - qsatt) - lstarn*dqsat - 1 + dlstarp*qp
            dT = -fff/dfff
            niter = niter + 1
            T1 = T1 + dT
        end

        qsatt = qsatt + dqsat * dT
        qn  = max(0.0, qt-qsatt)

    else
        qn = 0.0
    end
    
    T = T1#= - fac1*qp/(1+fac2*qp)
    if (T1 >= T0n)
        T = T1 - fac_cond*qp
    elseif (T1 <= T00n)
        T = T1 - fac_sub*qp
    end=#

    qp = max(0.0, qp)

    ωn = max(0,min(1,an*T-bn))
    ωp = max(0,min(1,ap*T-bp))
    ωg = max(0,min(1,ag*T-T00g*ag))

    qc = max(0.0,ωn*qn)
    qi = max(0.0,(1-ωn)*qn)

    qr = max(0.0,ωp*qp)
    qs = max(0.0,(1-ωp)*(1-ωg)*qp)
    qg = max(0.0,(1-ωp)*ωg*qp)
    P = moistPressure(PhysConst; ρ = ρ, Temp = T, qv = qt-qn)
    
    return T,P,qn,qc,qi,qr,qs,qg,qsatt
end


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
