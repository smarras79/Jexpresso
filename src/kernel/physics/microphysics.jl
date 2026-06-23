@kernel function do_micro_physics_gpu_3D!(uaux, qe, Tabs, qn, qi, qc, qr, qs, qg, Pr, Ps, Pg, S_micro, PhysConst, MicroConst, lpert, neq, npoin, z, adjusted, Pm)
    ip = @index(Global, Linear)

    T = eltype(uaux)
    if (lpert)
        @inbounds uaux[ip,6] = max(zero(T), (uaux[ip,6] + qe[ip,6])) - qe[ip,6]
        @inbounds uaux[ip,7] = max(zero(T), (uaux[ip,7] + qe[ip,7])) - qe[ip,7]
    else
        @inbounds uaux[ip,6] = max(zero(T), uaux[ip,6])
        @inbounds uaux[ip,7] = max(zero(T), uaux[ip,7])
    end
    uip = @view(uaux[ip,1:neq])
    qeip = @view(qe[ip,1:neq+1])

    @inbounds adjusted[ip,:] .= saturation_adjustment_sam_microphysics_gpu(uip,qeip,z[ip],MicroConst,PhysConst,lpert)
    @inbounds Tabs[ip] = adjusted[ip,1]
    @inbounds uaux[ip,end] = adjusted[ip,2]
    @inbounds qn[ip] = adjusted[ip,3]
    @inbounds qc[ip] = adjusted[ip,4]
    @inbounds qi[ip] = adjusted[ip,5]
    @inbounds qr[ip] = adjusted[ip,6]
    @inbounds qs[ip] = adjusted[ip,7]
    @inbounds qg[ip] = adjusted[ip,8]
    @inbounds qsatt = adjusted[ip,9]

    @inbounds Pm[ip,:] .= compute_Pm_gpu(uip, qeip, qr[ip], qs[ip], qg[ip], MicroConst,lpert)

    @inbounds Pr[ip] = Pm[ip,1]
    @inbounds Ps[ip] = Pm[ip,2]
    @inbounds Pg[ip] = Pm[ip,3]
    S = compute_dqpdt_sam_micro_gpu(uip, qeip, Tabs[ip], qn[ip], qc[ip], qi[ip], qr[ip], qs[ip], qg[ip], qsatt, MicroConst, PhysConst, lpert) 
    @inbounds S_micro[ip] = S
end

function do_micro_physics!(Tabs, qn, qc, qi, qr, qs, qg, Pr, Ps, Pg, S_micro, qsatt, npoin, uaux, z, qe, NSD, ::PERT)

     PhysConst = PhysicalConst{Float64}()
     MicroConst = MicrophysicalConst{Float64}()

    #get pressure and perform saturation adjustment
    saturation_adjustment_sam_microphysics!(uaux, qe, Tabs, qn, qi, qc, qr, qs, qg, qsatt, z, npoin, MicroConst, PhysConst, NSD, true)
    
    compute_Pm!(uaux, qe, qr, qs, qg, Pr, Ps, Pg, npoin, MicroConst, true)

    ### compute dqp/dt from microphysical processes
    compute_dqpdt_sam_micro!(uaux, qe, Tabs, qn, qc, qi, qr, qs, qg, qsatt, S_micro, npoin, MicroConst, PhysConst, NSD, true)

    ### find vertical derivatives for sedimentation processes.

    return nothing
end

function do_micro_physics!(Tabs, qn, qc, qi, qr, qs, qg, Pr, Ps, Pg, S_micro, qsatt, npoin, uaux, z, qe, NSD, ::TOTAL)

     PhysConst = PhysicalConst{Float64}()
     MicroConst = MicrophysicalConst{Float64}()

    #get pressure and perform saturation adjustment
    saturation_adjustment_sam_microphysics!(uaux, qe, Tabs, qn, qi, qc, qr, qs, qg, qsatt, z, npoin, MicroConst,PhysConst, NSD, false)

    compute_Pm!(uaux, qe, qr, qs, qg, Pr, Ps, Pg, npoin, MicroConst, false)

    ### compute dqp/dt from microphysical processes

    compute_dqpdt_sam_micro!(uaux, qe, Tabs, qn, qc, qi, qr, qs, qg, qsatt, S_micro, npoin, MicroConst, PhysConst, NSD, false)


    ### find vertical derivatives for sedimentation processes.

    return nothing
end

function compute_precipitation_derivatives!(dqpdt, dqtdt, dhldt, Pr, Ps, Pg, Tabs, qi, ρ, ρe, nelem, ngl, connijk, H, metrics, ω, dψ, SD::NSD_3D, ::TOTAL)

    MicroConst = MicrophysicalConst{Float64}()
    dqpdt .= 0.0
    dqtdt .= 0.0
    dhldt .= 0.0
    Lc     = MicroConst.Lc
    Ls     = MicroConst.Ls
    Lf     = MicroConst.Lf
    T0n    = MicroConst.T0n
    T00n   = MicroConst.T00n

    for e=1:nelem
        # do precipitation
        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip         = connijk[e,i,j,k]
                    H[i,j,k,1] = Pr[ip] + Ps[ip] + Pg[ip]
                end
            end
        end
        compute_vertical_derivative_q!(dqpdt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ,SD)
        
        # do precipitation effects on hl

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip         = connijk[e,i,j,k]
                    H[i,j,k,1] = Lc*Pr[ip] + Ls*(Ps[ip] + Pg[ip])
                end
            end
        end
        
        compute_vertical_derivative_q!(dhldt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ, SD)
        
        #do cloud ice sedimentation

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip = connijk[e,i,j,k]
                    H[i,j,k,1] = qi[ip] * ρ[ip] * 0.4 
                end
            end
        end
        
        compute_vertical_derivative_q!(dqtdt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ, SD)
        
        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip              = connijk[e,i,j,k]
                    T               = Tabs[ip]
                    ωn              = max(0,min(1,(T-T00n)/(T0n - T00n)))
                    dhldt[e,i,j,k] += (Lc + ωn*Lf)*dqtdt[e,i,j,k]
                end
            end
        end
    end

end

function compute_precipitation_derivatives!(dqpdt, dqtdt, dhldt, Pr, Ps, Pg, Tabs, qi, ρ, ρe, nelem, ngl, connijk, H, metrics, ω, dψ, SD::NSD_3D, ::PERT)

    MicroConst = MicrophysicalConst{Float64}()
    dqpdt .= 0.0
    dqtdt .= 0.0
    dhldt .= 0.0
    Lc = MicroConst.Lc
    Ls = MicroConst.Ls
    Lf = MicroConst.Lf
    T0n = MicroConst.T0n
    T00n = MicroConst.T00n
    Je = metrics.Je
    dξdz = metrics.dξdz
    dηdz = metrics.dηdz
    dζdz = metrics.dζdz
    for e=1:nelem
        # do precipitation
        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip         = connijk[e,i,j,k]
                    H[i,j,k,1] = Pr[ip] + Ps[ip] + Pg[ip]
                end
            end
        end
        compute_vertical_derivative_q!(dqpdt, H, e, ngl, Je, dξdz, dηdz, metrics.dζdz, ω, dψ, SD)

        # do precipitation effects on hl

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip         = connijk[e,i,j,k]
                    H[i,j,k,1] = Lc*Pr[ip] + Ls*(Ps[ip] + Pg[ip])
                end
            end
        end
        
        compute_vertical_derivative_q!(dhldt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ, SD)
        
        #do cloud ice sedimentation

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip         = connijk[e,i,j,k]
                    H[i,j,k,1] = qi[ip] * (ρ[ip] +ρe[ip])* 0.4
                end
            end
        end

        compute_vertical_derivative_q!(dqtdt, H, e, ngl, metrics.Je, metrics.dξdz, metrics.dηdz, metrics.dζdz,ω,dψ, SD)

        for i=1:ngl
            for j=1:ngl
                for k=1:ngl
                    ip              = connijk[e,i,j,k]
                    T               = Tabs[ip]
                    ωn              = max(0,min(1,(T-T00n)/(T0n - T00n)))
                    dhldt[e,i,j,k] += (Lc + ωn*Lf)*dqtdt[e,i,j,k]
                end
            end
        end
    end

end

function compute_precipitation_derivatives!(dqpdt, dqtdt, dhldt, Pr, Ps, Pg, Tabs, qi, ρ, ρe, nelem, ngl, connijk, H, metrics, ω, dψ, SD::NSD_2D, ::TOTAL)

    MicroConst = MicrophysicalConst{Float64}()
    dqpdt .= 0.0
    dqtdt .= 0.0
    dhldt .= 0.0
    Lc = MicroConst.Lc
    Ls = MicroConst.Ls
    Lf = MicroConst.Lf
    T0n = MicroConst.T0n
    T00n = MicroConst.T00n
    for e=1:nelem
        # do precipitation
        for i=1:ngl
            for j=1:ngl
                ip         = connijk[e,i,j]
                H[i,j,1] = Pr[ip] + Ps[ip] + Pg[ip]
            end
        end
        compute_vertical_derivative_q!(dqpdt, H, e, ngl, metrics.Je, metrics.dξdy, metrics.dηdy, ω, dψ, SD)

        # do precipitation effects on hl

        for i=1:ngl
            for j=1:ngl
                ip         = connijk[e,i,j]
                H[i,j,1] = Lc*Pr[ip] + Ls*(Ps[ip] + Pg[ip])
            end
        end

        compute_vertical_derivative_q!(dhldt, H, e, ngl, metrics.Je, metrics.dξdy, metrics.dηdy, ω,dψ, SD)

        #do cloud ice sedimentation

        for i=1:ngl
            for j=1:ngl
                ip         = connijk[e,i,j]
                H[i,j,1] = qi[ip] * (ρ[ip])* 0.4
            end
        end
        compute_vertical_derivative_q!(dqtdt, H, e, ngl, metrics.Je, metrics.dξdy, metrics.dηdy, ω,dψ, SD)

        for i=1:ngl
            for j=1:ngl
                ip              = connijk[e,i,j]
                T               = Tabs[ip]
                ωn              = max(0,min(1,(T-T00n)/(T0n - T00n)))
                dhldt[e,i,j] += (Lc + ωn*Lf)*dqtdt[e,i,j]
            end
        end
    end

end


function compute_precipitation_derivatives!(dqpdt, dqtdt, dhldt, Pr, Ps, Pg, Tabs, qi, ρ, ρe, nelem, ngl, connijk, H, metrics, ω, dψ, SD::NSD_2D, ::PERT)

    MicroConst = MicrophysicalConst{Float64}()
    dqpdt .= 0.0
    dqtdt .= 0.0
    dhldt .= 0.0
    Lc = MicroConst.Lc
    Ls = MicroConst.Ls
    Lf = MicroConst.Lf
    T0n = MicroConst.T0n
    T00n = MicroConst.T00n
    for e=1:nelem
        # do precipitation
        for i=1:ngl
            for j=1:ngl
                @inbounds ip::Int64         = connijk[e,i,j]
                @inbounds H[i,j,1] = Pr[ip] + Ps[ip] + Pg[ip]
            end
        end
        compute_vertical_derivative_q!(dqpdt, H, e, ngl, metrics.Je, metrics.dξdy, metrics.dηdy, ω, dψ, SD)

        # do precipitation effects on hl

        for i=1:ngl
            for j=1:ngl
                @inbounds ip         = connijk[e,i,j]
                @inbounds H[i,j,1] = Lc*Pr[ip] + Ls*(Ps[ip] + Pg[ip])
            end
        end

        compute_vertical_derivative_q!(dhldt, H, e, ngl, metrics.Je, metrics.dξdy, metrics.dηdy, ω,dψ, SD)

        #do cloud ice sedimentation

        for i=1:ngl
            for j=1:ngl
                @inbounds ip         = connijk[e,i,j]
                @inbounds H[i,j,1] = qi[ip] * (ρ[ip] +ρe[ip])* 0.4
            end
        end

        compute_vertical_derivative_q!(dqtdt, H, e, ngl, metrics.Je, metrics.dξdy, metrics.dηdy, ω,dψ, SD)

        for i=1:ngl
            for j=1:ngl
                @inbounds ip              = connijk[e,i,j]
                @inbounds T               = Tabs[ip]
                ωn              = max(0,min(1,(T-T00n)/(T0n - T00n)))
                @inbounds dhldt[e,i,j] += (Lc + ωn*Lf)*dqtdt[e,i,j]
            end
        end
    end

end


function precipitation_flux_gpu(u,qe,MicroConst,lpert,Pr,Ps,Pg,qi)
    T= eltype(u)
    if (lpert)
        ρ = u[1] + qe[1]
    else
        ρ = u[1]
    end
    Lc = MicroConst.Lc
    Ls = MicroConst.Ls

    return T(0.0), T(Lc*Pr +Ls*(Ps+Pg)), T(qi*ρ*T(0.4)), T(Pr + Ps + Pg) 
end

function add_micro_precip_sources!(S, q, qe, S_micro::Float64, qn::Float64, flux_lw::Float64, flux_sw::Float64, PhysConst, ::NSD_3D, ::TOTAL)
    
    @inbounds ρ        = q[1]
    @inbounds qt       = q[6]/ρ
    @inbounds qp       = q[7]/ρ
    qv::Float64        = qt - qn
    @inbounds ρqv_pert = ρ*qv - qe[6]
     
    @inbounds S[4] += PhysConst.g*(0.61*ρqv_pert -ρ*(qn+qp))
    @inbounds S[5] += ρ*(flux_lw - flux_sw)
    @inbounds S[6] += -ρ*S_micro
    @inbounds S[7] += ρ*S_micro

end

function add_micro_precip_sources!(S, q, qe, S_micro::Float64, qn::Float64, flux_lw::Float64, flux_sw::Float64, PhysConst, ::NSD_3D, ::PERT)
    
    @inbounds ρ        = q[1]+qe[1]
    @inbounds qt       = (q[6] + qe[6])/ρ
    @inbounds qv_ref   = qe[6]/qe[1]
    @inbounds qp       = (q[7] + qe[7]) /ρ
    qv::Float64        = (qt - qn) #- qe[6]/qe[1]
    @inbounds ρqv_pert = ρ*qv - qe[6]
    
    @inbounds S[4] += PhysConst.g*(0.61*ρqv_pert -ρ*(qn+qp))# should we ignore condensates in the hydrostatic balance if they're not included in the pressure term?
    @inbounds S[5] += ρ*(flux_lw - flux_sw)
    @inbounds S[6] += -ρ*S_micro
    @inbounds S[7] += ρ*S_micro

end

function add_micro_precip_sources!(S, q, qe, S_micro::Float64, qn::Float64, flux_lw::Float64, flux_sw::Float64, PhysConst, ::NSD_2D, ::TOTAL)   

    
    @inbounds ρ        = q[1]
    @inbounds qt       = q[5]/ρ
    @inbounds qp       = q[6]/ρ
    qv::Float64        = qt - qn
    @inbounds ρqv_pert = ρ*qv - qe[5]

    @inbounds S[3] += PhysConst.g*(0.61*ρqv_pert -ρ*(qn+qp)) 
    @inbounds S[4] += ρ*(flux_lw - flux_sw)
    @inbounds S[5] += -ρ*S_micro
    @inbounds S[6] += ρ*S_micro

end 
    
function add_micro_precip_sources!(S, q, qe, S_micro::Float64, qn::Float64, flux_lw::Float64, flux_sw::Float64, PhysConst, ::NSD_2D, ::PERT)
    
    
    @inbounds ρ        = q[1]+qe[1]
    @inbounds ρqt      = (q[5] + qe[5])
    qt::Float64        = ρqt/ρ
    @inbounds qv_ref   = qe[5]/qe[1]
    @inbounds ρqp      = (q[6] + qe[6])
    qp::Float64        = ρqp/ρ
    qv::Float64        = (qt - qn) #- qe[6]/qe[1]
    @inbounds ρqv_pert = ρ*qv - qe[5]

    @inbounds S[3] += PhysConst.g*(0.61*ρqv_pert -ρ*(qn+qp))# should we ignore condensates in the hydrostatic balance if they're not included in the pressure term?
    @inbounds S[4] += ρ*(flux_lw - flux_sw)
    @inbounds S[5] += -ρ*S_micro
    @inbounds S[6] += ρ*S_micro

end

function precipitation_source_gpu(u,qe,lpert,qn,S_micro,PhysConst,MicroConst)
    T = eltype(u)
    if (lpert)
        ρ   = u[1] + qe[1]
        qt  = (u[6] + qe[6])/ρ
        qp  = (u[7] + qe[7])/ρ
        ρqv = ρ*(qt - qn) - qe[6]
    else
        ρ   = u[1]
        qt  = u[6]/ρ
        qp  = u[7]/ρ
        ρqv = ρ*(qt - qn)
    end
    
    return T(-PhysConst.g*(T(0.608)*ρqv-ρ*(qn+qp))), T(0.0), T(-ρ*S_micro), T(ρ*S_micro)
end

function compute_dqpdt_sam_micro!(uaux, qe, Tabs, qn, qc, qi, qr, qs, qg, qsatt, S_micro, npoin, MicroConst, PhysConst, NSD, lpert)
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
    for ip=1:npoin
        if (lpert)
            if (NSD == NSD_3D())
                ρ = uaux[ip,1] + qe[ip,1]
                qt = (uaux[ip,6] + qe[ip,6])/ρ
                qp = (uaux[ip,7] + qe[ip,7])/ρ
            else
                ρ = uaux[ip,1] + qe[ip,1]
                qt = (uaux[ip,5] + qe[ip,5])/ρ
                qp = (uaux[ip,6] + qe[ip,6])/ρ
            end
        else
            if (NSD == NSD_3D())
                ρ = uaux[ip,1]
                qt = uaux[ip,6]/uaux[ip,1]
                qp = uaux[ip,7]/uaux[ip,1]
            else
                ρ = uaux[ip,1]
                qt = uaux[ip,5]/uaux[ip,1]
                qp = uaux[ip,6]/uaux[ip,1]
            end

        end
        T = Tabs[ip]
        e_satw = esatw(T)*100
        e_sati = esati(T)*100

        ### Collection of condensates
        Ar_c = π/4*a_rain*N0_rain*Er_c*γ3br*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_rain*N0_rain))^((3+b_rain)/4)
        Ar_i = exp(0.025*(T-273.16))*π/4*a_rain*N0_rain*Er_i*γ3br*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_rain*N0_rain))^((3+b_rain)/4)
        As_c = π/4*a_snow*N0_snow*Es_c*γ3bs*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_snow*N0_snow))^((3+b_snow)/4)
        As_i = exp(0.025*(T-273.16))*π/4*a_snow*N0_snow*Es_i*γ3bs*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_snow*N0_snow))^((3+b_snow)/4)
        Ag_c = π/4*a_graupel*N0_graupel*Eg_c*γ3bs*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_graupel*N0_graupel))^((3+b_graupel)/4)
        Ag_i = exp(0.025*(T-273.16))*π/4*a_graupel*N0_graupel*Eg_c*γ3bs*(ρ0/ρ)^(0.5)*(ρ/(π*ρ_graupel*N0_graupel))^((3+b_graupel)/4)

        dqrdt = (Ar_c*qc[ip] + Ar_i*qi[ip])*qr[ip]^((3+b_rain)/4)
        dqsdt = (As_c*qc[ip] + As_i*qi[ip])*qs[ip]^((3+b_snow)/4)
        dqgdt = (Ag_c*qc[ip] + Ag_i*qi[ip])*qg[ip]^((3+b_graupel)/4)

        ### Autoconversion

        Auto = max(0.0, α*(qc[ip] - qc0))

        ### Aggregation

        Aggr = max(0.0, β*exp(0.025*(T-273.16))*(qi[ip]-qi0))

        ### Evaporation

        S = max((qt-qn[ip]),0.0)/qsatt[ip]
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

        Evap_r = min(0.0,(2*π*C_rain*N0_rain)/(ρ*(A_rain+B_rain))*(A_er*sqrt(qr[ip]) + B_er*qr[ip]^((5+b_rain)/8))*(S-1))
        Evap_s = min(0.0,(2*π*C_snow*N0_snow)/(ρ*(A_snow_graupel+B_snow_graupel))*(A_es*sqrt(qs[ip]) + B_es*qs[ip]^((5+b_snow)/8))*(S-1))
        Evap_g = min(0.0,(2*π*C_graupel*N0_graupel)/(ρ*(A_snow_graupel+B_snow_graupel))*(A_eg*sqrt(qg[ip]) + B_eg*qg[ip]^((5+b_graupel)/8))*(S-1))

        ### cloud ice sedimentation is handled as part of the precipitation package, this routine ignores it

        S_micro[ip] = clamp(dqrdt + dqsdt + dqgdt + Auto + Aggr + Evap_r + Evap_s + Evap_g, -qp, qn[ip])
    end
end


function compute_dqpdt_sam_micro_gpu(u,qe,T,qn,qc,qi,qr,qs,qg,qsatt,MicroConst,PhysConst,lpert)
    FT = eltype(u)
    if (lpert)
            ρ = u[1] + qe[1]
            qt = (u[6] + qe[6])/ρ
            qp = (u[7] + qe[7])/ρ
            P = u[end]
    else
            ρ = u[1]
            qt = u[6]/ρ
            qp = u[7]/ρ
            P = u[end]
    end

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
    e_satw = esatw(T)*FT(100)
    e_sati = esati(T)*FT(100)

    ### Collection of condensates
    Ar_c = FT(FT(π)/FT(4))*a_rain*N0_rain*Er_c*γ3br*(ρ0/ρ)^(FT(0.5))*FT((ρ/(FT(π)*ρ_rain*N0_rain)))^((FT(3)+FT(b_rain)/FT(4)))
    Ar_i = exp(FT(0.025)*(T-FT(273.16)))*FT(FT(π)/FT(4))*a_rain*N0_rain*Er_i*γ3br*(FT(ρ0/ρ))^(FT(0.5))*FT((ρ/(FT(π)*ρ_rain*N0_rain)))^(FT((FT(3)+b_rain)/FT(4)))
    As_c = FT(FT(π)/FT(4))*a_snow*N0_snow*Es_c*γ3bs*(FT(ρ0/ρ))^(FT(0.5))*(FT(ρ/(FT(π)*ρ_snow*N0_snow)))^(FT((FT(3)+b_snow)/FT(4)))
    As_i = exp(FT(0.025)*(T-FT(273.16)))*FT(FT(π)/FT(4))*a_snow*N0_snow*Es_i*γ3bs*(FT(ρ0/ρ))^(FT(0.5))*FT((ρ/(FT(π)*ρ_snow*N0_snow)))^(FT((FT(3)+b_snow)/FT(4)))
    Ag_c = FT(FT(π)/FT(4))*a_graupel*N0_graupel*Eg_c*γ3bs*(FT(ρ0/ρ))^(FT(0.5))*FT((ρ/(FT(π)*ρ_graupel*N0_graupel)))^(FT((FT(3)+b_graupel)/FT(4)))
    Ag_i = exp(FT(0.025)*(T-FT(273.16)))*FT(FT(π)/FT(4))*a_graupel*N0_graupel*Eg_c*γ3bs*(FT(ρ0/ρ))^(FT(0.5))*(FT(ρ/(FT(π)*ρ_graupel*N0_graupel)))^(FT((FT(3)+b_graupel)/FT(4)))

    dqrdt = (Ar_c*qc + Ar_i*qi)*qr^(FT((FT(3)+b_rain)/FT(4)))
    dqsdt = (As_c*qc + As_i*qi)*qs^(FT((FT(3)+b_snow)/FT(4)))
    dqgdt = (Ag_c*qc + Ag_i*qi)*qg^(FT((FT(3)+b_graupel)/FT(4)))

    ### Autoconversion

    Auto = FT(max(FT(0.0), α*(qc - qc0)))

    ### Aggregation

    Aggr = FT(max(FT(0.0), β*exp(FT(0.025)*(T-FT(273.16)))*(qi-qi0)))

    ### Evaporation

    S = FT((qt-qn)/qsatt)
    A_rain = FT((Lc/(Ka*T)))*FT(((Lc/(Rvap*T))) - FT(1))
    A_snow_graupel = FT((Ls/(Ka*T)))*FT(((Ls/(Rvap*T))) - FT(1))
    B_rain = FT(Rvap*Rair/(Da*e_satw))
    B_snow_graupel = FT(Rvap*Rair/(Da*e_sati))
    A_er = a_fr*FT((ρ/(FT(π)*ρ_rain*N0_rain)))^(FT(0.5))
    A_es = a_fs*FT((ρ/(FT(π)*ρ_snow*N0_snow)))^(FT(0.5))
    A_eg = a_fg*FT((ρ/(FT(π)*ρ_graupel*N0_graupel)))^(FT(0.5))
    B_er = b_fr*FT((ρ*a_rain/μ))^(FT(0.5))*γ5br*FT((ρ0/ρ))^(FT(0.25))*FT((ρ/(FT(π)*ρ_rain*N0_rain)))^(FT((FT(5)+b_rain)/FT(8)))
    B_es = b_fs*FT((ρ*a_snow/μ))^(FT(0.5))*γ5bs*FT((ρ0/ρ))^(FT(0.25))*FT((ρ/(FT(π)*ρ_snow*N0_snow)))^(FT((FT(5)+b_snow)/FT(8)))
    B_eg = b_fg*FT((ρ*a_graupel/μ))^(FT(0.5))*γ5bg*FT((ρ0/ρ))^(FT(0.25))*FT((ρ/(FT(π)*ρ_graupel*N0_graupel)))^(FT((FT(5)+b_graupel)/FT(8)))

    Evap_r = FT((FT(2)*FT(π)*C_rain*N0_rain)/(ρ*(A_rain+B_rain)))*(A_er*sqrt(qr) + B_er*qr^(FT((FT(5)+b_rain)/FT(8))))*(S-FT(1))
    Evap_s = FT((FT(2)*FT(π)*C_snow*N0_snow)/(ρ*(A_snow_graupel+B_snow_graupel)))*(A_es*sqrt(qs) + B_es*qs^(FT((FT(5)+b_snow)/FT(8))))*(S-FT(1))
    Evap_g = FT((FT(2)*FT(π)*C_graupel*N0_graupel)/(ρ*(A_snow_graupel+B_snow_graupel)))*(A_eg*sqrt(qg) + B_eg*qg^(FT((FT(5)+b_graupel)/FT(8))))*(S-FT(1))

    ### cloud ice sedimentation is handled as part of the precipitation package, this routine ignores it

    dqpdt = FT(clamp(dqrdt + dqsdt + dqgdt + Auto + Aggr + Evap_r + Evap_s + Evap_g,-qp,qn))

    #if (qn > 0.001)
        #@info dqrdt, dqsdt, dqgdt, Auto, Aggr, Evap_r, Evap_s, Evap_g, qn, qc, qi, qp, qr, qs, qg, dqpdt
    #end

    return FT(dqpdt)
end

function compute_Pm!(uaux, qe, qr, qs, qg, Pr, Ps, Pg, npoin, MicroConst, lpert)
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
    for ip=1:npoin
        if (lpert)
            ρ = uaux[ip,1] + qe[ip,1] 
        else
            ρ = uaux[ip,1]
        end
        #@info ρ,qr[ip],ρ0
        Pr[ip] = (a_rain * γ4br)/6 * (π * ρ_rain * N0_rain)^(-b_rain/4)*(ρ0/ρ)^(0.5)*(ρ*qr[ip])^(1+b_rain/4)
        Ps[ip] = (a_snow * γ4bs)/6 * (π * ρ_snow * N0_snow)^(-b_snow/4)*(ρ0/ρ)^(0.5)*(ρ*qs[ip])^(1+b_snow/4)
        Pg[ip] = (a_graupel * γ4bg)/6 * (π * ρ_graupel * N0_graupel)^(-b_graupel/4)*(ρ0/ρ)^(0.5)*(ρ*qg[ip])^(1+b_graupel/4)
    end
end

function compute_Pm_gpu(u, qe, qr, qs, qg, MicroConst, lpert)
   
    FT = eltype(u)
    if (lpert)
        ρ = u[1] + qe[1]
    else
        ρ = u[1]
    end
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

    Pr = (a_rain * γ4br)/FT(6) * (FT(π) * ρ_rain * N0_rain)^(FT(-b_rain/FT(4)))*(ρ0/ρ)^(FT(0.5))*(ρ*qr)^(FT(1)+FT(b_rain/FT(4)))
    Ps = (a_snow * γ4bs)/FT(6) * (FT(π) * ρ_snow * N0_snow)^(FT(-b_snow/FT(4)))*(ρ0/ρ)^(FT(0.5))*(ρ*qs)^(FT(1)+FT(b_snow/FT(4)))
    Pg = (a_graupel * γ4bg)/FT(6) * (FT(π) * ρ_graupel * N0_graupel)^(FT(-b_graupel/FT(4)))*(ρ0/ρ)^(FT(0.5))*(ρ*qg)^(FT(1)+FT(b_graupel/FT(4)))
    #if (Pr > 0)
    #    @info Pr, Ps, Pg
    #end
    return FT(Pr), FT(Ps), FT(Pg)
end

function saturation_adjustment_sam_microphysics!(uaux, qe, Tabs, qn, qi, qc, qr, qs, qg, qsatt, z, npoin, MicroConst, PhysConst, NSD, lpert)

    T00n = MicroConst.T00n #Temperature threshold for cloud water
    T0n  = MicroConst.T0n #Temperature threshold for ice
    T00p = MicroConst.T00p #Temperature threshold for rain     
    T0p  = MicroConst.T0p #Temperature threshold for snow/graupel 
    T00g = MicroConst.T00g #Temperature threshold for graupel
    T0g  = MicroConst.T0g #Temperature threshold for graupel
    
    Lc = MicroConst.Lc
    Lf = MicroConst.Lf
    Ls = MicroConst.Ls
    
    g  = PhysConst.g
    cp = PhysConst.cp
    
    fac_cond = Lc/cp
    fac_fus  = Lf/cp
    fac_sub  = Ls/cp

    for ip=1:npoin
        if (lpert)
            if (NSD == NSD_3D())
                ρ  = uaux[ip,1] + qe[ip,1]
                hl = (uaux[ip,5] + qe[ip,5])/ρ
                qt = (uaux[ip,6] + qe[ip,6])/ρ
                qp = (uaux[ip,7] + qe[ip,7])/ρ
                qt = max(0.0, qt)
                qp = max(0.0, qp)
            
                uaux[ip,6] = qt*ρ - qe[ip,6]
                uaux[ip,7] = qp*ρ - qe[ip,7]
            else

                ρ  = uaux[ip,1] + qe[ip,1]
                hl = (uaux[ip,4] + qe[ip,4])/ρ
                qt = (uaux[ip,5] + qe[ip,5])/ρ
                qp = (uaux[ip,6] + qe[ip,6])/ρ
                qt = max(0.0, qt)
                qp = max(0.0, qp)

                uaux[ip,5] = qt*ρ - qe[ip,5]
                uaux[ip,6] = qp*ρ - qe[ip,6]

            end
        
        else
            if (NSD == NSD_3D())
                ρ  = uaux[ip,1]
                hl = uaux[ip,5]/ρ
                qt = max(0.0,uaux[ip,6])/ρ
                qp = max(0.0,uaux[ip,7])/ρ
            
                uaux[ip,6] = max(0.0,uaux[ip,6])
                uaux[ip,7] = max(0.0,uaux[ip,7])
            else
                ρ  = uaux[ip,1]
                hl = uaux[ip,4]/ρ
                qt = max(0.0,uaux[ip,5])/ρ
                qp = max(0.0,uaux[ip,6])/ρ

                uaux[ip,5] = max(0.0,uaux[ip,5])
                uaux[ip,6] = max(0.0,uaux[ip,6])
            end
        end

        # find equilibrium temperature from saturation adjustment
        # initial guess for sensible temperature and pressure assumes no condensates/all vapor
                
        T =  (hl - g*z[ip])/cp
        an   = 1/(T0n - T00n)
        bn   = T00n * an
        ap   = 1/(T0p - T00p)
        bp   = T00p*ap
        fac1 = fac_cond+(1+bp)*fac_fus
        fac2 = fac_fus*ap
        ag   = 1/(T0g - T00g)
        ωp   = max(0,min(1,ap*T-bp))
        T1   = T + (fac_cond + (1-ωp)*fac_fus)*qp #+ fac1*qp/(1+fac2*qp)
        Tv   = T1*(1 + 0.61*qt - qp)
        P    = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)
        if (T1 >= T0p)
            
            T1 = T + fac_cond*qp
            Tv = T1*(1 + 0.61*qt - qp)
            P  = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)
        
        elseif (T1 <= T00p)
            
            T1 = T + fac_sub*qp
            Tv = T1*(1 + 0.61*qt - qp)
            P  = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)
        
        else
            
            ωp = max(0,min(1,ap*T1-bp))
            T1 = T + (fac_cond + (1-ωp)*fac_fus)*qp
            Tv = T1*(1 + 0.61*qt - qp)
            P  = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)
        
        end
        if (T1 >= T0n)

            qsatt[ip] = max(0.0,qsatw(T1, P/100))
        
        elseif (T1 <= T00n)

            qsatt[ip] = max(0.0,qsati(T1, P/100))
        
        else
        
            ωn        = max(0,min(1,an*T1-bn))
            qsatt[ip] = max(0.0,ωn*qsatw(T1,P/100)+(1-ωn)*qsati(T1,P/100))
        
        end
        
        Tabs[ip] = T1

        if (qt > qsatt[ip])
            Tv = T1*(1 + 0.61*min(qt,qsatt[ip]) - qp - max(0,qt-qsatt[ip]))
            P  = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)
            
            niter = 0
            dT    = 100
            dqsat = 0.0
            
            while (abs(dT) > 0.0001 && niter < 50)

                if (T1 >= T0n)

                    ωn        = 1
                    lstarn    = fac_cond
                    dlstarn   = 0
                    qsatt[ip] = max(qsatw(T1,P/100),0.0)
                    dqsat     = dtqsatw(T1,P/100)
                    dωn       = 0.0
                
                elseif (T1 <= T00n)

                    ωn        = 0
                    lstarn    = fac_sub
                    dlstarn   = 0
                    qsatt[ip] = max(qsati(T1,P/100),0.0)
                    dqsat     = dtqsati(T1,P/100)
                    dωn       = 0.0
                
                else

                    ωn        = max(0,min(1,an*T1-bn))
                    dωn       = an
                    lstarn    = fac_cond+(1-ωn)*fac_fus
                    dlstarn   = -dωn*fac_fus#dωn*fac_cond - dωn * fac_fus 
                    qsatt[ip] = max(0.0,ωn*qsatw(T1,P/100) + (1-ωn)*qsati(T1,P/100))
                    dqsat     = ωn*dtqsati(T1,P/100) + (1-ωn)*dtqsati(T1,P/100) + dωn * qsatw(T1,P/100) - dωn * qsati(T1,P/100)
                
                end

                if (T1 >= T0p)

                    ωp      = 1
                    lstarp  = fac_cond
                    dlstarp = 0

                elseif (T1 <= T00p)

                    ωp      = 0
                    lstarp  = fac_sub
                    dlstarp = 0

                else

                    ωp      = max(0,min(1,ap*T1-bp))
                    lstarp  = fac_cond + (1-ωp)*fac_fus
                    dlstarp = -ap*fac_fus#ap*fac_cond - ap*fac_fus#ap*fac_fus

                end

                fff   = T - T1 + lstarn*(qt - qsatt[ip]) + lstarp*qp
                dfff  = dlstarn*(qt - qsatt[ip]) - lstarn*dqsat - 1 + dlstarp*qp
                dT    = -fff/dfff
                niter = niter + 1
                T1    = T1 + dT
                Tv    = T1*(1 + 0.61*min(qt,qsatt[ip]) - qp - max(0,qt-qsatt[ip]))
                P     = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)
            
            end
            
            qsatt[ip] = qsatt[ip] + dT*dqsat
            qn[ip]    = max(0.0, qt-qsatt[ip])

        else
        
            qn[ip] = 0.0
        
        end

        Tabs[ip] = T1

        ωn = max(0,min(1,an*Tabs[ip]-bn))
        ωp = max(0,min(1,ap*Tabs[ip]-bp))
        ωg = max(0,min(1,ag*Tabs[ip]-T00g*ag))

        qc[ip] = max(0.0,ωn*qn[ip])
        qi[ip] = max(0.0,(1-ωn)*qn[ip])

        qr[ip] = max(0.0,ωp*qp)
        qs[ip] = max(0.0,(1-ωp)*(1-ωg)*qp)
        qg[ip] = max(0.0,(1-ωp)*ωg*qp)
        
        uaux[ip,end] = moistPressure(PhysConst; ρ = ρ, Tv = Tv, qv = qt-qn[ip])
    end
end


function saturation_adjustment_sam_microphysics_gpu(u,qe,z,MicroConst,PhysConst,lpert)

    FT = eltype(u)
    if (lpert)
        
        ρ = u[1] + qe[1]
        hl = (u[5] + qe[5])/ρ
        qt = (u[6] + qe[6])/ρ
        qp = (u[7] + qe[7])/ρ
    
    else
    
        ρ = u[1]
        hl = u[5]/ρ
        qt = u[6]/ρ
        qp = u[7]/ρ
    
    end

    T00n = MicroConst.T00n #Temperature threshold for cloud water
    T0n  = MicroConst.T0n #Temperature threshold for ice
    T00p = MicroConst.T00p #Temperature threshold for rain     
    T0p  = MicroConst.T0p #Temperature threshold for snow/graupel 
    T00g = MicroConst.T00g #Temperature threshold for graupel
    T0g  = MicroConst.T0g #Temperature threshold for graupel

    Lc = MicroConst.Lc
    Lf = MicroConst.Lf
    Ls = MicroConst.Ls
    
    g  =PhysConst.g
    cp = PhysConst.cp
    
    fac_cond = Lc/cp
    fac_fus  = Lf/cp
    fac_sub  = Ls/cp
    # find equilibrium temperature from saturation adjustment
    # initial guess for sensible temperature and pressure assumes no condensates/all vapor
    
    T =  (hl - g*z)/cp
    qt = FT(max(FT(0.0),qt))
    qp = FT(max(FT(0.0),qp))
    #P = moistPressure(PhysConst; ρ=ρ, Temp=T, qv = qt) 
    
    an   = FT(1/(T0n - T00n))
    bn   = T00n * an
    ap   = FT(1/(T0p - T00p))
    bp   = T00p*ap
    fac1 = fac_cond+(FT(1)+bp)*fac_fus
    fac2 = fac_fus*ap
    ag   = FT(FT(1)/(T0g - T00g))
    T1   = T + fac1*qp/(FT(1)+fac2*qp)
   
    Tv   = T1*(1 + 0.61*qt - qp)
    P    = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)

    if (T1 >= T0p)

        T1 = T + fac_cond*qp
        Tv = T1*(1 + 0.61*qt - qp)
        P  = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)

    elseif (T1 <= T00p)

        T1 = T + fac_sub*qp
        Tv = T1*(1 + 0.61*qt - qp)
        P  = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)

    else

        ωp = max(0,min(1,ap*T1-bp))
        T1 = T + (fac_cond + (1-ωp)*fac_fus)*qp
        Tv = T1*(1 + 0.61*qt - qp)
        P  = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)

    end

    if (T1 >= T0n)

        qsatt = FT(max(0.0,qsatw(T1, P/100)))

    elseif (T1 <= T00n)

        qsatt = FT(max(0.0,qsati(T1, P/100)))

    else

        ωn    = FT(max(0,min(1,an*T1-bn)))
        qsatt = FT(max(0.0,ωn*qsatw(T1,P/100)+(1-ωn)*qsati(T1,P/100)))

    end

    if (qt > qsatt)

        niter = Int32(0)
        dT    = FT(100)
        dqsat = FT(0.0)
        
        while (abs(dT) > FT(0.001) && niter < FT(50))

            if (T1 >= T0n)

                ωn      = 1
                lstarn  = fac_cond
                dlstarn = FT(0)
                qsatt   = FT(qsatw(T1,P/100))
                dqsat   = FT(dtqsatw(T1,P/100))

            elseif (T1 <= T00n)

                ωn      = FT(0)
                lstarn  = fac_sub
                dlstarn = FT(0)
                qsatt   = FT(qsati(T1,P/100))
                dqsat   = FT(dtqsati(T1,P/100))

            else

                ωn      = FT(max(FT(0),FT(min(FT(1),an*T1-bn))))
                lstarn  = fac_cond+(FT(1)-ωn)*fac_fus
                dlstarn = an*fac_fus
                qsatt   = FT(max(FT(0.0),ωn*FT(qsatw(T1,P/100)) + (FT(1)-ωn)*FT(qsati(T1,P/100))))
                dqsat   = ωn*FT(dtqsati(T1,P/100)) + (FT(1)-ωn)*FT(dtqsati(T1,P/100)) + dωn * qsatw(T1,P/100) - dωn * qsati(T1,P/100)
            end

            if (T1 >= T0p)

                ωp      = FT(1)
                lstarp  = fac_cond
                dlstarp = FT(0)

            elseif (T1 <= T00p)

                ωp      = FT(0)
                lstarp  = fac_sub
                dlstarp = FT(0)

            else

                ωp      = FT(max(FT(0),FT(min(FT(1),ap*T1-bp))))
                lstarp  = fac_cond + (FT(1)-ωp)*fac_fus
                dlstarp = ap*fac_fus

            end

            fff   = T - T1 + lstarn*(qt - qsatt) + lstarp*qp
            dfff  = dlstarn*(qt - qsatt) - lstarn*dqsat - FT(1) + dlstarp*qp
            dT    = -fff/dfff
            niter = niter + Int32(1)
            T1    = T1 + dT
            Tv    = T1*(FT(1) + FT(0.61)*min(qt,qsatt) - qp - max(FT(0),qt-qsatt))
            P     = moistPressure(PhysConst; ρ=ρ, Tv=Tv, qv = qt)
        end

        qsatt = qsatt + dqsat * dT
        qn    = max(FT(0.0), qt-qsatt)

    else

        qn = FT(0.0)
    end

    T = T1

    qp = FT(max(FT(0.0), qp))
    ωn = FT(max(FT(0),FT(min(FT(1),an*T-bn))))
    ωp = FT(max(FT(0),FT(min(FT(1),ap*T-bp))))
    ωg = FT(max(FT(0),FT(min(FT(1),ag*T-T00g*ag))))

    qc = FT(max(FT(0.0),ωn*qn))
    qi = FT(max(FT(0.0),(FT(1)-ωn)*qn))

    qr = FT(max(FT(0.0),ωp*qp))
    qs = FT(max(FT(0.0),(FT(1)-ωp)*(FT(1)-ωg)*qp))
    qg = FT(max(FT(0.0),(FT(1)-ωp)*ωg*qp))
    P = moistPressure(PhysConst; ρ = ρ, Temp = T, qv = qt-qn)

    return FT(T),FT(P),FT(qn),FT(qc),FT(qi),FT(qr),FT(qs),FT(qg),FT(qsatt)
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
                      @view(q[:,1]),          #ρ
                      @view(q[:,mesh.nsd+2]), #θ
                      @view(q[:,mesh.nsd+2]), #qv
                      @view(q[:,mesh.nsd+3]), #qc
                      @view(q[:,mesh.nsd+4]), #qr
                      qref,
                      @view(mesh.coords[:,end]),
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
    qrr, pressure, temp,                         = 0.0, 0.0, 0.0
    es, qvs, dz, pii                             = 0.0, 0.0, 0.0, 0.0
    f5, dtfall, rdz, product                     = 0.0, 0.0, 0.0, 0.0
    max_heating, max_condense, max_rain, maxqrp  = 0.0, 0.0, 0.0, -100.0
    vtmax, ernmax, crmax, factorn                = 0.0, 0.0, 0.0, 0.0
    time_sediment                                = dt
    qcr, factorr, ppt                            = 0.0, 0.0, 0.0
    max_cr_sedimentation                         = 0.75
    
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
