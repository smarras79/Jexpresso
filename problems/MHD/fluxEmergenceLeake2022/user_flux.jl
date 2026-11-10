#---------------------------------------------------------------------------------
# 3D flux emergence into a coronal point dipole, the parameter study of Leake,
# Linton & Antiochos, "The role of reconnection in the onset of solar eruptions",
# ApJ 934:10 (2022). Ideal GLM-MHD flux in three dimensions, with γ = 5/3 (the same
# flux as problems/MHD/fluxEmergenceLeake2014).
#
# Units (Sec. 2.2): L₀ = 170 km, B₀ = 0.13 T, g₀ = g_sun = 274 m s⁻², hence
# P₀ = B₀²/μ₀ = 1.345e4 Pa, ρ₀ = P₀/(L₀g₀) = 2.887e-4 kg m⁻³, v₀ = √(L₀g₀) =
# 6.825 km s⁻¹, t₀ = 24.909 s, T₀ = m_p L₀ g₀/k_B = 5637 K, J₀ = B₀/(μ₀L₀) =
# 0.609 A m⁻², E₀ = ρ₀v₀²L₀³ = 6.607e19 J. In these units g = 1, the magnetic pressure is
# ½|B|² (μ₀ is absorbed, as in Heaviside-Lorentz units), and the ideal gas law
# P = ρ k_B T/(m_f m_p) reads p = ρT/m_f; the paper's atmosphere has m_f = 1 in these
# units (initialize.jl, README).
#
# State (ρ, ρu, ρv, E, ρw, Bx, By, Bz, ψ): slot 4 is the total energy and slot 5
# the vertical momentum (z is up), the layout of every MHD case of this tree.
#---------------------------------------------------------------------------------
if !@isdefined(γ_mhd)
    const γ_mhd = 5.0/3.0
end
if abs(γ_mhd - 5.0/3.0) > 1e-12
    error(" problems/MHD/fluxEmergenceLeake2022: γ_mhd = $(γ_mhd) is already defined by another MHD case loaded in this Julia session (this case needs γ = 5/3). Restart Julia before running this case.")
end
if !@isdefined(c_h_mhd)
    const c_h_mhd = Ref{Float64}(1.0)
end

@inline function pressure_mhd(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz, ψ)
    γm1 = γ_mhd - 1.0
    ke  = 0.5*(ρu*ρu + ρv*ρv + ρw*ρw)/ρ
    me  = 0.5*(Bx*Bx + By*By + Bz*Bz)
    return γm1*(ρE - ke - me - 0.5*ψ*ψ)
end

function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=9, ip=1)

    ρ  = q[1]; ρu = q[2]; ρv = q[3]; ρE = q[4]; ρw = q[5]
    Bx = q[6]; By = q[7]; Bz = q[8]; ψ  = q[9]
    u = ρu/ρ; v = ρv/ρ; w = ρw/ρ
    c_h = c_h_mhd[]

    p    = pressure_mhd(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz, ψ)
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
    # y
    G[1] = ρv
    G[2] = ρu*v - By*Bx
    G[3] = ρv*v + ptot - By*By
    G[4] = v*enfl - By*vdB + c_h*ψ*By
    G[5] = ρw*v - By*Bz
    G[6] = v*Bx - u*By
    G[7] = c_h*ψ
    G[8] = v*Bz - w*By
    G[9] = c_h*By
    # z
    H[1] = ρw
    H[2] = ρu*w - Bz*Bx
    H[3] = ρv*w - Bz*By
    H[4] = w*enfl - Bz*vdB + c_h*ψ*Bz
    H[5] = ρw*w + ptot - Bz*Bz
    H[6] = w*Bx - u*Bz
    H[7] = w*By - v*Bz
    H[8] = c_h*ψ
    H[9] = c_h*Bz
end

function user_flux!(F, G, H, q, qe, mesh::St_mesh, ::CL, ::PERT; neqs=9, ip=1)
    error(" problems/MHD/fluxEmergenceLeake2022: PERT() solution variables are not supported.")
end
