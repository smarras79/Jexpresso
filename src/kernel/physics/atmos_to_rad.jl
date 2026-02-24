" Convert atmospheric state to absorption and scattering coefficients"

function atmos_to_rad(atmos_data,npoin)

    PhysConst = PhysicalConst{Float64}()

    σ = zeros(npoin)
    κ = zeros(npoin)
    for ip = 1:npoin

        T = atmos_data.t_lay[ip]
        P = atmos_data.p_lay[ip]
        qv = atmos_data.vmr_h2o[ip]
        q_liq = atmos_data.q_liq[ip]
        q_ice = atmos_data.q_ice[ip]
        ρ = atmos_data.rho[ip]

        β = 1.66
        m_l = 1.2
        m_i = 0.9
        
        c = β * m_l
        d = β * m_i
        z1 = (PhysConst.Rair/PhysConst.g)*T*log(PhysConst.pref/20)
        b = 2.5/z1

        a = 0.0001

        κ[ip] = ρ*(a + b*qv + c*q_liq + d*q_ice)
        σ[ip] = 0.1*κ[ip]
    end
    @info " done computing extinction and scattering here are extrema", extrema(κ), extrema(σ)
    return κ, σ
end
