#---------------------------------------------------------------------------------
# 2.5D visco-resistive GLM-MHD: emergence of a twisted flux tube from the
# convection zone into a field-free, quadrupolar or dipolar corona,
#
#   J. E. Leake, M. G. Linton & S. K. Antiochos, "Tests of dynamical flux
#   emergence as a mechanism for coronal mass ejection initiation",
#   ApJ 722:550 (2010), https://doi.org/10.1088/0004-637X/722/1/550,
#   Erratum: ApJ 741:125 (2011), https://doi.org/10.1088/0004-637X/741/2/125.
#
# The paper writes its equations in Gaussian units with dimensional values (Mm,
# s, G, g cm⁻³). The code units are photospheric, as in the Lare runs of the
# same authors (Arber, Haynes & Leake 2007, ApJ 666:541: "normalised by ...
# photospheric values"; the erratum quotes L0 = 150 km, v0 = 6.5 km/s):
#
#   L0 = 150 km             length
#   T0 = T_ph = 5700 K      temperature           (T_ph = 1)
#   v0 = (R T_ph)^½         velocity, 6.467 km/s  (R = k_B/(μ m_p), μ = 1.125)
#   t0 = L0/v0 = 23.19 s    time (the "23 s" of the paper's Eq. 6)
#   ρ0 = ρ_ph, p0 = ρ0 v0² = p_ph                 (ρ_ph = p_ph = 1)
#   B0 = (4π p_ph)^½ = 1030.6 G                   magnetic pressure ½|B|²
#
# so that the gas law is p = ρT, gravity is g = g_⊙ L0/v0² = 0.98273, the
# photospheric pressure scale height is T_ph/g = 1.0176 L0 = 152.6 km. The
# mean molecular weight μ and the photospheric pressure p_ph are not printed in
# the paper; both are fitted to the initial profiles of its Fig. 3(a) (see
# README.md, "Units and the reference atmosphere").
#
# State (ρ, ρu, ρv, E, ρw, Bx, By, Bz, ψ): the layout of every MHD case of this
# tree (slot 4 the total energy, slot 5 the out-of-plane momentum). Jexpresso's
# (x, y) is the paper's (x, y), y up; z is the paper's ignorable direction, so
# w and Bz are its shear flow and shear field.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Physical constants and the paper's parameters (cgs, Gaussian)
#---------------------------------------------------------------------------------
if !@isdefined(LK10_L0)
    const LK10_kB    = 1.380649e-16      # erg K⁻¹
    const LK10_mp    = 1.67262192e-24    # g
    const LK10_clight= 2.99792458e10     # cm s⁻¹
    const LK10_L0    = 1.5e7             # cm, the erratum's normalizing length (150 km)
    const LK10_Tph   = 5700.0            # K, photosphere/chromosphere (Sec. 2.2.1)
    const LK10_Tcor  = 8.6e5             # K, corona
    const LK10_ytr   = 3.75e8            # cm, height of the transition region
    const LK10_wtr   = 0.75e8            # cm, its width
    const LK10_gsun  = 2.74e4            # cm s⁻², |g| (Leake & Linton 2013 state 274 m s⁻² for these runs)
    const LK10_MU    = 1.125             # mean molecular weight, fitted to Fig. 3(a)
    const LK10_PPH   = 10.0^4.927        # dyn cm⁻², p(y = 0), fitted to Fig. 3(a)
    # derived code units
    const LK10_V0    = sqrt(LK10_kB*LK10_Tph/(LK10_MU*LK10_mp))       # cm s⁻¹
    const LK10_T0    = LK10_L0/LK10_V0                                 # s
    const LK10_RHO0  = LK10_PPH/LK10_V0^2                              # g cm⁻³ (ρ_ph)
    const LK10_BU    = sqrt(4π*LK10_PPH)                               # G
    const LK10_G     = LK10_gsun*LK10_L0/LK10_V0^2                     # code gravity
    const LK10_TCOR  = LK10_Tcor/LK10_Tph                              # code T_cor
    const LK10_YTR   = LK10_ytr/LK10_L0                                # 25
    const LK10_WTR   = LK10_wtr/LK10_L0                                # 5
    # Viscosity ν = 2.6e4 g cm⁻¹ s⁻¹ (Sec. 2.1), DYNAMIC: ν/(ρ_ph L0 v0) in code
    # units. Used here as a KINEMATIC viscosity at that value everywhere (README,
    # "Deviations", item 2).
    const LK10_NU    = 2.6e4/(LK10_RHO0*LK10_L0*LK10_V0)              # 0.01326
    # Anomalous resistivity (Eqs. 7-8): η = η₀ max(0, |j|/j_crit − 1) with
    # η₀ = 1.6e-11 s and j_crit = 2.07e7 statA cm⁻². As a magnetic diffusivity
    # D = ηc²/4π, and |j| = c|∇×B|/4π, in code units:
    const LK10_ETA0  = 1.6e-11*LK10_clight^2/(4π)/(LK10_L0*LK10_V0)   # 1.18e-4
    const LK10_JCRIT = 4π*2.07e7/LK10_clight*LK10_L0/LK10_BU          # 126.3 (|∇×B| in B0/L0)
end

#---------------------------------------------------------------------------------
# Ideal GLM-MHD flux, γ = 5/3 (Sec. 2.1: "polytropic index γ = 5/3")
#---------------------------------------------------------------------------------
# γ_mhd, c_h_mhd and pressure_mhd are shared by name with the other MHD cases of
# this tree (one equation set). A Julia session that already ran a γ ≠ 5/3 case
# keeps that value, hence the check.
if !@isdefined(γ_mhd)
    const γ_mhd = 5.0/3.0
end
if abs(γ_mhd - 5.0/3.0) > 1e-12
    error(" problems/MHD/fluxEmergenceLeake2010: γ_mhd = $(γ_mhd) is already defined by another MHD case loaded in this Julia session (this case needs γ = 5/3). Restart Julia before running this case.")
end

# Hyperbolic divergence-cleaning speed c_h, set by initialize.jl.
if !@isdefined(c_h_mhd)
    const c_h_mhd = Ref{Float64}(1.0)
end

# Pressure floor used ONLY inside the flux evaluation (as in
# problems/MHD/fluxEmergenceSon2025): 1% of the smallest pressure of the initial
# atmosphere (p = 1.0e-7 at y = 600), never active in a healthy run.
if !@isdefined(LK10_PFLOOR)
    const LK10_PFLOOR = Ref{Float64}(1.0e-9)
end

@inline function pressure_mhd(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz, ψ)
    γm1 = γ_mhd - 1.0
    ke  = 0.5*(ρu*ρu + ρv*ρv + ρw*ρw)/ρ
    me  = 0.5*(Bx*Bx + By*By + Bz*Bz)
    return γm1*(ρE - ke - me - 0.5*ψ*ψ)
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=9, ip=1)

    ρ  = q[1]; ρu = q[2]; ρv = q[3]; ρE = q[4]; ρw = q[5]
    Bx = q[6]; By = q[7]; Bz = q[8]; ψ  = q[9]
    u = ρu/ρ; v = ρv/ρ; w = ρw/ρ
    c_h = c_h_mhd[]

    p    = max(pressure_mhd(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz, ψ), LK10_PFLOOR[])
    magp = 0.5*(Bx*Bx + By*By + Bz*Bz)
    ptot = p + magp
    vdB  = u*Bx + v*By + w*Bz
    enfl = 0.5*(ρu*u + ρv*v + ρw*w) + γ_mhd*p/(γ_mhd - 1.0) + 2.0*magp   # E + p + ½|B|² − ½ψ²

    # x
    F[1] = ρu
    F[2] = ρu*u + ptot - Bx*Bx
    F[3] = ρv*u - Bx*By
    F[4] = u*enfl - Bx*vdB + c_h*ψ*Bx
    F[5] = ρw*u - Bx*Bz
    F[6] = c_h*ψ
    F[7] = u*By - v*Bx
    F[8] = u*Bz - w*Bx
    F[9] = c_h*Bx
    # y (vertical)
    G[1] = ρv
    G[2] = ρu*v - By*Bx
    G[3] = ρv*v + ptot - By*By
    G[4] = v*enfl - By*vdB + c_h*ψ*By
    G[5] = ρw*v - By*Bz
    G[6] = v*Bx - u*By
    G[7] = c_h*ψ
    G[8] = v*Bz - w*By
    G[9] = c_h*By
end
