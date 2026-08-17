function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=5, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρe = q[5]

    u  = ρu/ρ
    v  = ρv/ρ
    w  = ρw/ρ

    γ         = PhysConst.γ
    γm1       = γ - 1.0
    velomagsq = u*u + v*v + w*w
    ke        = 0.5*ρ*velomagsq
    Pressure  = γm1*(ρe - ke)

    # Total-energy enthalpy flux: H_flux = u*(ρe + p) = u*(ke + γ*p/γm1)
    F[1] = ρu
    F[2] = ρu*u + Pressure
    F[3] = ρv*u
    F[4] = ρw*u
    F[5] = u*(ke + γ*Pressure/γm1)

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Pressure
    G[4] = ρw*v
    G[5] = v*(ke + γ*Pressure/γm1)

    H[1] = ρw
    H[2] = ρu*w
    H[3] = ρv*w
    H[4] = ρw*w + Pressure
    H[5] = w*(ke + γ*Pressure/γm1)
end

@inline convert_transformed_to_primitive(u, ::GradientVariablesPrimitive) = u
@inline convert_derivative_to_primitive(u, gradient, ::GradientVariablesPrimitive) = gradient

@inline function convert_transformed_to_primitive(w, ::GradientVariablesEntropy)
    PhysConst = PhysicalConst{Float64}()
    u   = entropy2cons(w)
    rho = u[1]
    v1  = u[2]/rho; v2 = u[3]/rho; v3 = u[4]/rho
    p   = PhysConst.γm1 * (u[5] - 0.5*rho*(v1*v1 + v2*v2 + v3*v3))
    T   = p / (rho * PhysConst.Rair)
    return SVector(rho, v1, v2, v3, T)
end

@inline function convert_derivative_to_primitive(w, grad_w, ::GradientVariablesEntropy)
    # w2=ρv1/p, w3=ρv2/p, w4=ρv3/p, w5=-ρ/p ⇒ p/ρ = -1/w5.
    #   ∂vi = (p/ρ)(∂w_{i+1} + vi ∂w5),  ∂T = (1/Rair)(p/ρ)² ∂w5
    PhysConst = PhysicalConst{Float64}()
    invR = 1.0 / PhysConst.Rair
    w5   = w[5]
    pρ   = -1.0 / w5
    v1   = -w[2]/w5; v2 = -w[3]/w5; v3 = -w[4]/w5
    return SVector(grad_w[1],
                   pρ * (grad_w[2] + v1*grad_w[5]),
                   pρ * (grad_w[3] + v2*grad_w[5]),
                   pρ * (grad_w[4] + v3*grad_w[5]),
                   invR * pρ * pρ * grad_w[5])
end

@inline function cons2entropy(u)
    PhysConst = PhysicalConst{Float64}()
    γ      = PhysConst.γ
    invγm1 = 1.0 / PhysConst.γm1
    rho, rho_v1, rho_v2, rho_v3, rho_e = u[1], u[2], u[3], u[4], u[5]
    v1 = rho_v1/rho; v2 = rho_v2/rho; v3 = rho_v3/rho
    v_square = v1*v1 + v2*v2 + v3*v3
    p = PhysConst.γm1 * (rho_e - 0.5*rho*v_square)
    s = log(p) - γ*log(rho)
    rho_p = rho/p
    w1 = (γ - s)*invγm1 - 0.5*rho_p*v_square
    w2 = rho_p*v1
    w3 = rho_p*v2
    w4 = rho_p*v3
    w5 = -rho_p
    return SVector(w1, w2, w3, w4, w5)
end

@inline function entropy2cons(w)
    PhysConst = PhysicalConst{Float64}()
    γ      = PhysConst.γ
    γm1    = PhysConst.γm1
    invγm1 = 1.0 / γm1
    V1, V2, V3, V4, V5 = w[1]*γm1, w[2]*γm1, w[3]*γm1, w[4]*γm1, w[5]*γm1
    s = γ - V1 + (V2*V2 + V3*V3 + V4*V4)/(2.0*V5)
    rho_iota = (γm1 / (-V5)^γ)^invγm1 * exp(-s*invγm1)
    rho    = -rho_iota*V5
    rho_v1 =  rho_iota*V2
    rho_v2 =  rho_iota*V3
    rho_v3 =  rho_iota*V4
    rho_e  =  rho_iota*(1.0 - (V2*V2 + V3*V3 + V4*V4)/(2.0*V5))
    return SVector(rho, rho_v1, rho_v2, rho_v3, rho_e)
end

# Effective dynamic viscosity μ for the NS parabolic flux, dispatched on the viscosity model:
#   AV   → constant molecular μ = visc_coeffieq[2]      (DNS reference)
#   SMAG → Smagorinsky (μ_mol+μ_turb)·visc_coeffieq[2]  (LES); visc_coeffieq[2] is the SGS
#          intensity multiplier swept over {1.0,0.5,0.25}.
@inline ns_effective_mu(::AV, vc, rho, g11,g22,g33,g12,g21,g13,g31,g23,g32, PhysConst, Δ2, inputs) =
    vc[2]

# LES: μ_eff = μ_mol (fixed, = Re=1600 value via inputs[:mu_molecular]) + intensity·μ_turb,
# with the Smagorinsky eddy viscosity μ_turb = ρ·C_s²·Δ²·|S| computed directly so the molecular
# floor is decoupled from the SGS-intensity sweep (vc[2] = JEXP_MU scales only μ_turb).
@inline function ns_effective_mu(::SMAG, vc, rho, g11,g22,g33,g12,g21,g13,g31,g23,g32, PhysConst, Δ2, inputs)
    Cs  = get(inputs, :C_smag, PhysConst.C_s)::Float64
    S12 = 0.5*(g12+g21); S13 = 0.5*(g13+g31); S23 = 0.5*(g23+g32)
    SijSij = g11*g11 + g22*g22 + g33*g33 + 2.0*(S12*S12 + S13*S13 + S23*S23)
    Smag = sqrt(2.0*SijSij)
    mu_turb = rho * Cs*Cs * Δ2 * Smag
    mu_mol  = inputs[:mu_molecular]::Float64
    return mu_mol + vc[2]*mu_turb
end

# 3D Navier–Stokes viscous (parabolic) flux. `u` are the transformed gradient variables
# (primitive ρ,v1,v2,v3,T or entropy w1..w5); `gradients`=(grad_x,grad_y,grad_z) of them.
function flux_parabolic(u, gradients, orientation::Integer, visc_coeffieq, inputs, delta2, gradvars, VT)
    rho, v1, v2, v3, _ = convert_transformed_to_primitive(u, gradvars)
    _, dv1dx, dv2dx, dv3dx, dTdx = convert_derivative_to_primitive(u, gradients[1], gradvars)
    _, dv1dy, dv2dy, dv3dy, dTdy = convert_derivative_to_primitive(u, gradients[2], gradvars)
    _, dv1dz, dv2dz, dv3dz, dTdz = convert_derivative_to_primitive(u, gradients[3], gradvars)

    # Newtonian deviatoric stress (τ_ii = 2∂vi - (2/3)∇·v); symmetric off-diagonals.
    div3   = dv1dx + dv2dy + dv3dz
    tau_11 = 2.0*dv1dx - (2.0/3.0)*div3
    tau_22 = 2.0*dv2dy - (2.0/3.0)*div3
    tau_33 = 2.0*dv3dz - (2.0/3.0)*div3
    tau_12 = dv1dy + dv2dx
    tau_13 = dv1dz + dv3dx
    tau_23 = dv2dz + dv3dy

    PhysConst = PhysicalConst{Float64}()
    gamma = PhysConst.γ
    Pr    = inputs[:Pr]::Float64

    mu = ns_effective_mu(VT, visc_coeffieq, rho,
                         dv1dx, dv2dy, dv3dz, dv1dy, dv2dx, dv1dz, dv3dx, dv2dz, dv3dy,
                         PhysConst, delta2, inputs)

    # Fourier heat flux, applied as (...+q)·mu ⇒ k = mu·cp/Pr, cp = γ·Rair/(γ-1).
    # kappa carries Rair because uprimitive[5] = T = p/(ρ·Rair) (Jexpresso R = Rair).
    kappa = gamma * PhysConst.Rair / ((gamma - 1) * Pr)
    q1 = kappa*dTdx; q2 = kappa*dTdy; q3 = kappa*dTdz

    if orientation == 1
        return SVector(0.0, tau_11*mu, tau_12*mu, tau_13*mu,
                       (v1*tau_11 + v2*tau_12 + v3*tau_13 + q1)*mu)
    elseif orientation == 2
        return SVector(0.0, tau_12*mu, tau_22*mu, tau_23*mu,
                       (v1*tau_12 + v2*tau_22 + v3*tau_23 + q2)*mu)
    else
        return SVector(0.0, tau_13*mu, tau_23*mu, tau_33*mu,
                       (v1*tau_13 + v2*tau_23 + v3*tau_33 + q3)*mu)
    end
end

@inline ln_mean(x::Real, y::Real) = ln_mean(promote(x, y)...)
@inline function ln_mean(x::RealT, y::RealT) where {RealT <: Real}
    eps_f2 = convert(RealT, 1.0e-4)
    f2 = (x*(x - 2*y) + y*y) / (x*(x + 2*y) + y*y)
    f2 < eps_f2 ? (x + y)/@evalpoly(f2, 2, convert(RealT,2/3), convert(RealT,2/5), convert(RealT,2/7)) :
                  (y - x)/log(y/x)
end

@inline inv_ln_mean(x::Real, y::Real) = inv_ln_mean(promote(x, y)...)
@inline function inv_ln_mean(x::RealT, y::RealT) where {RealT <: Real}
    eps_f2 = convert(RealT, 1.0e-4)
    f2 = (x*(x - 2*y) + y*y) / (x*(x + 2*y) + y*y)
    f2 < eps_f2 ? @evalpoly(f2, 2, convert(RealT,2/3), convert(RealT,2/5), convert(RealT,2/7))/(x + y) :
                  log(y/x)/(y - x)
end

@inline function user_fluxaux!(aux, SD::NSD_3D, q, ::TOTAL,
                               ::Union{ranocha,kennedy_gruber,shima,chandrashekar})
    PhysConst = PhysicalConst{Float64}()
    rho = q[1]
    ru = q[2]
    rv = q[3]
    rw = q[4]
    re = q[5]
    u = ru/rho
    v = rv/rho
    w = rw/rho
    p = PhysConst.γm1*(re - 0.5*(ru*u + rv*v + rw*w))
    aux[1]=rho
    aux[2]=u
    aux[3]=v
    aux[4]=w
    aux[5]=p
    aux[6]=re
end

@inline function flux_turbo(u_ll, u_rr, ::ranocha)
    γ = PhysicalConst{Float64}().γ
    rho_ll,v1_ll,v2_ll,v3_ll,p_ll,_ = u_ll
    rho_rr,v1_rr,v2_rr,v3_rr,p_rr,_ = u_rr
    rho_mean = ln_mean(rho_ll, rho_rr)
    inv_rho_p_mean = p_ll*p_rr*inv_ln_mean(rho_ll*p_rr, rho_rr*p_ll)
    v1_avg=0.5*(v1_ll+v1_rr)
    v2_avg=0.5*(v2_ll+v2_rr)
    v3_avg=0.5*(v3_ll+v3_rr)
    p_avg=0.5*(p_ll+p_rr)
    vsq=0.5*(v1_ll*v1_rr+v2_ll*v2_rr+v3_ll*v3_rr)
    ig=1.0/(γ-1.0)
    f1=rho_mean*v1_avg
    g1=rho_mean*v2_avg
    h1=rho_mean*v3_avg
    eint = vsq + inv_rho_p_mean*ig
    F=SVector(f1, f1*v1_avg+p_avg, f1*v2_avg, f1*v3_avg, f1*eint+0.5*(p_ll*v1_rr+p_rr*v1_ll))
    G=SVector(g1, g1*v1_avg, g1*v2_avg+p_avg, g1*v3_avg, g1*eint+0.5*(p_ll*v2_rr+p_rr*v2_ll))
    H=SVector(h1, h1*v1_avg, h1*v2_avg, h1*v3_avg+p_avg, h1*eint+0.5*(p_ll*v3_rr+p_rr*v3_ll))
    return F, G, H
end

@inline function flux_turbo(u_ll, u_rr, ::kennedy_gruber)
    rho_ll,v1_ll,v2_ll,v3_ll,p_ll,re_ll = u_ll
    rho_rr,v1_rr,v2_rr,v3_rr,p_rr,re_rr = u_rr
    rho_avg=0.5*(rho_ll+rho_rr)
    v1_avg=0.5*(v1_ll+v1_rr)
    v2_avg=0.5*(v2_ll+v2_rr)
    v3_avg=0.5*(v3_ll+v3_rr)
    p_avg=0.5*(p_ll+p_rr)
    e_avg=0.5*(re_ll/rho_ll + re_rr/rho_rr)
    hh=rho_avg*e_avg+p_avg
    f1=rho_avg*v1_avg
    g1=rho_avg*v2_avg
    h1=rho_avg*v3_avg
    F=SVector(f1, f1*v1_avg+p_avg, f1*v2_avg, f1*v3_avg, hh*v1_avg)
    G=SVector(g1, g1*v1_avg, g1*v2_avg+p_avg, g1*v3_avg, hh*v2_avg)
    H=SVector(h1, h1*v1_avg, h1*v2_avg, h1*v3_avg+p_avg, hh*v3_avg)
    return F, G, H
end

@inline function flux_turbo(u_ll, u_rr, ::shima)
    γ = PhysicalConst{Float64}().γ
    rho_ll,v1_ll,v2_ll,v3_ll,p_ll,_ = u_ll
    rho_rr,v1_rr,v2_rr,v3_rr,p_rr,_ = u_rr
    rho_avg=0.5*(rho_ll+rho_rr)
    v1_avg=0.5*(v1_ll+v1_rr)
    v2_avg=0.5*(v2_ll+v2_rr)
    v3_avg=0.5*(v3_ll+v3_rr)
    p_avg=0.5*(p_ll+p_rr)
    kin=0.5*(v1_ll*v1_rr+v2_ll*v2_rr+v3_ll*v3_rr)
    ig=1.0/(γ-1.0)
    f1=rho_avg*v1_avg
    g1=rho_avg*v2_avg
    h1=rho_avg*v3_avg
    F=SVector(f1, f1*v1_avg+p_avg, f1*v2_avg, f1*v3_avg, p_avg*v1_avg*ig+f1*kin+0.5*(p_ll*v1_rr+p_rr*v1_ll))
    G=SVector(g1, g1*v1_avg, g1*v2_avg+p_avg, g1*v3_avg, p_avg*v2_avg*ig+g1*kin+0.5*(p_ll*v2_rr+p_rr*v2_ll))
    H=SVector(h1, h1*v1_avg, h1*v2_avg, h1*v3_avg+p_avg, p_avg*v3_avg*ig+h1*kin+0.5*(p_ll*v3_rr+p_rr*v3_ll))
    return F, G, H
end

@inline function flux_turbo(u_ll, u_rr, ::chandrashekar)
    γ = PhysicalConst{Float64}().γ
    rho_ll,v1_ll,v2_ll,v3_ll,p_ll,_ = u_ll
    rho_rr,v1_rr,v2_rr,v3_rr,p_rr,_ = u_rr
    beta_ll=0.5*rho_ll/p_ll; beta_rr=0.5*rho_rr/p_rr
    vsq=0.5*(v1_ll^2+v2_ll^2+v3_ll^2) + 0.5*(v1_rr^2+v2_rr^2+v3_rr^2)
    rho_avg=0.5*(rho_ll+rho_rr); rho_mean=ln_mean(rho_ll,rho_rr)
    beta_mean=ln_mean(beta_ll,beta_rr); beta_avg=0.5*(beta_ll+beta_rr)
    v1_avg=0.5*(v1_ll+v1_rr); v2_avg=0.5*(v2_ll+v2_rr); v3_avg=0.5*(v3_ll+v3_rr)
    p_mean=0.5*rho_avg/beta_avg; ef=0.5*(1.0/((γ-1.0)*beta_mean) - vsq)
    f1=rho_mean*v1_avg; f2=f1*v1_avg+p_mean; f3=f1*v2_avg; f4=f1*v3_avg
    g1=rho_mean*v2_avg; g2=g1*v1_avg; g3=g1*v2_avg+p_mean; g4=g1*v3_avg
    h1=rho_mean*v3_avg; h2=h1*v1_avg; h3=h1*v2_avg; h4=h1*v3_avg+p_mean
    F=SVector(f1,f2,f3,f4, f1*ef + f2*v1_avg + f3*v2_avg + f4*v3_avg)
    G=SVector(g1,g2,g3,g4, g1*ef + g2*v1_avg + g3*v2_avg + g4*v3_avg)
    H=SVector(h1,h2,h3,h4, h1*ef + h2*v1_avg + h3*v2_avg + h4*v3_avg)
    return F, G, H
end
