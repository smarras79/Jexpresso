#---------------------------------------------------------------------------------
# 2.5D analogue of the flux-emergence study of Leake, Linton & Antiochos, ApJ 934:10
# (2022): ideal GLM-MHD flux in the plane, all three components of v and B, γ = 5/3.
#
# Units as problems/MHD/fluxEmergenceLeake2022 (L₀ = 170 km, B₀ = 0.13 T, g₀ = 274
# m s⁻², hence v₀ = 6.82 km s⁻¹, t₀ = 24.9 s): g = 1, magnetic pressure ½|B|², and
# p = ρT/m_f with m_f = 1 (README).
#
# Coordinates: the code's 2D plane is (x, y) with y UP (the paper's z); the rope axis,
# the paper's y, is the invariant direction. State (ρ, ρu, ρv, E, ρw, Bx, By, Bz, ψ):
# u, Bx horizontal, v, By vertical, w, Bz along the rope (out of the plane).
#---------------------------------------------------------------------------------
if !@isdefined(γ_mhd)
    const γ_mhd = 5.0/3.0
end
if abs(γ_mhd - 5.0/3.0) > 1e-12
    error(" problems/MHD/fluxEmergenceLeake2022_2D: γ_mhd = $(γ_mhd) is already defined by another MHD case loaded in this Julia session (this case needs γ = 5/3). Restart Julia before running this case.")
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

function user_flux!(F, G, SD::NSD_2D,
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

    # x (horizontal)
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

function user_flux!(F, G, SD::NSD_2D, q, qe, mesh::St_mesh, ::CL, ::PERT; neqs=9, ip=1)
    error(" problems/MHD/fluxEmergenceLeake2022_2D: PERT() solution variables are not supported.")
end
