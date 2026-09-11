#---------------------------------------------------------------------------------
# 1D ideal MHD fluxes for the Brio-Wu shock tube (Brio & Wu 1988; Dao &
# Nazarov 2022, JSC 92:77, Sec. 5.2).
#
# Unknowns, in Jexpresso's MHD slot order minus the GLM field:
#
#     q = (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz)
#
# with Bx constant (∂x Bx = 0 is the whole of ∇·B = 0 in 1D, so no cleaning
# is needed; the Bx slot carries a zero flux and stays at its initial
# value). Magnetic pressure ½|B|², total pressure p_T = p + ½|B|², γ = 2.
#
#     F = ( ρu,
#           ρu² + p_T − Bx²,
#           ρuv − Bx By,
#           (ρE + p_T) u − Bx (u Bx + v By + w Bz),
#           ρuw − Bx Bz,
#           0,
#           By u − Bx v,
#           Bz u − Bx w )
#
# p = (γ − 1)(ρE − ½ρ|v|² − ½|B|²), floored at p_floor_mhd inside the flux
# only (the DynSGS regularization keeps it positive in a healthy run).
#---------------------------------------------------------------------------------
if !@isdefined(γ_mhd)
    const γ_mhd = 2.0
end
if abs(γ_mhd - 2.0) > 1e-12
    error(" problems/MHD/brioWu1d: γ_mhd = $(γ_mhd) is already defined by another MHD case loaded in this Julia session (this case needs γ = 2). Restart Julia before running this case.")
end
if !@isdefined(p_floor_mhd)
    const p_floor_mhd = 1.0e-9
end

@inline function pressure_mhd1d(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz)
    ke = 0.5*(ρu*ρu + ρv*ρv + ρw*ρw)/ρ
    me = 0.5*(Bx*Bx + By*By + Bz*Bz)
    return max((γ_mhd - 1.0)*(ρE - ke - me), p_floor_mhd)
end

function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=8, ip=1)

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρE = q[4]
    ρw = q[5]
    Bx = q[6]
    By = q[7]
    Bz = q[8]

    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ
    p  = pressure_mhd1d(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz)
    pT = p + 0.5*(Bx*Bx + By*By + Bz*Bz)
    uB = u*Bx + v*By + w*Bz

    F[1] = ρu
    F[2] = ρu*u + pT - Bx*Bx
    F[3] = ρv*u - Bx*By
    F[4] = (ρE + pT)*u - Bx*uB
    F[5] = ρw*u - Bx*Bz
    F[6] = 0.0
    F[7] = By*u - Bx*v
    F[8] = Bz*u - Bx*w
end
