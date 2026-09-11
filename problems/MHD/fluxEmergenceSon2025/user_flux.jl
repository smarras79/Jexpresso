#---------------------------------------------------------------------------------
# 2D ideal GLM-MHD equations for the flux-emergence problem of
#
#   D. Son, Y. Jang, T. Magara,
#   "A Comparative Analysis of High-resolution Shock-capturing Schemes for
#    Two-dimensional Magnetohydrodynamic Simulation of Flux Emergence in the
#    Solar Atmosphere", ApJS 277:46 (2025), Section 2.2, Eqs. (6)-(8).
#   https://doi.org/10.3847/1538-4365/adb617
#
# This is the SAME conservative GLM-MHD flux already used by
# problems/MHD/orszagTangBormanis2024 and problems/MHD/kelvinHelmholtzChan2022,
# with two differences dictated by the paper:
#
#   1. γ = 1.05 (Section 2.2: "taken to be 1.05 in this paper. This reduced
#      value of γ is chosen as it leads to a higher growth rate of the undular
#      mode of the magnetic buoyancy instability and mimics quasi-isothermal
#      heating with limited energy injection in the coronal temperature").
#   2. Heaviside-Lorentz units (the paper's choice, its Section 2.2): the
#      magnetic pressure is ½|B|², there are no 4π factors anywhere. The
#      paper's Eq. (2) is written in Gaussian units, B = [8πp/β]^½; in the
#      code units used here it reads B = [2p/β]^½, see initialize.jl.
#
# State (2D in space; the third components of velocity and magnetic field are
# carried for generality of the equation set and stay identically zero here):
#
#   q = (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)
#
# NOTICE on the ordering: the paper writes U = (ρ, ρV, B, e, ψ). Here the
# total energy ρE sits in slot 4 and the out-of-plane momentum ρw in slot 5
# because Jexpresso's shared 2D kernels (viscous τ·u augmentation, DynSGS
# residual floors, sound-speed/CFL diagnostic) assume that slot 4 of a 2D
# system carries the energy.
#
# Jexpresso's y is the paper's vertical coordinate z: gravity acts along -y
# (user_source.jl), "v" is the paper's V_z and "By" its B_z.
#
# Pressure:  p = (γ-1)·(ρE - ½ρ|v|² - ½|B|² - ½ψ²)          (paper Eq. 8)
#
# The GLM fluxes are Dedner's; the ψ damping (paper Eq. 7, last row) and the
# gravity sources (rows 3 and 8 of the paper's S) live in user_source.jl.
# The non-conservative Powell/Galilean-GLM term is not needed (no
# split-form / entropy-stable discretization here) and is omitted, as in the
# two companion MHD cases.
#---------------------------------------------------------------------------------

# Heat-capacity ratio of the paper. These constants are shared by name with
# the other MHD cases (they implement one and the same equation set), so a
# Julia session that has already run one of them keeps THEIR value — a
# γ = 5/3 leak into this γ = 1.05 case would be silent, hence the check.
if !@isdefined(γ_mhd)
    const γ_mhd = 1.05
end
if abs(γ_mhd - 1.05) > 1e-12
    error(" problems/MHD/fluxEmergenceSon2025: γ_mhd = $(γ_mhd) is already defined by another MHD case loaded in this Julia session (this case needs γ = 1.05). Restart Julia before running this case.")
end

# Nondimensional gravity g₀ = C_s²/(γ H₀) (paper Table 1): with H₀ = C_s = 1
# this is 1/γ. Points along -y.
if !@isdefined(g_mhd)
    const g_mhd = 1.0/γ_mhd
end

# Hyperbolic divergence-cleaning speed c_h. Filled by initialize.jl with the
# maximum wave speed |v| + c_f of the initial condition (paper Eq. 11) and
# kept constant throughout the run.
if !@isdefined(c_h_mhd)
    const c_h_mhd = Ref{Float64}(1.0)
end

# Pressure floor used ONLY inside the flux evaluation. In the corona the
# plasma β drops to ~10⁻⁴ (paper Fig. 6(e)) and, with γ - 1 = 0.05, the gas
# pressure is a 10⁻⁵ residue of the total energy: p = 0.05·(E - ½ρ|v|² -
# ½|B|²). A discretization error of the large terms can then push p through
# zero. The floor keeps the flux finite; it is 1% of the smallest initial
# pressure of the atmosphere (p ≈ 1.6e-7 at z = 35 H₀) and never active in a
# healthy run. Tunable from the REPL (Ref).
if !@isdefined(p_floor_mhd)
    const p_floor_mhd = Ref{Float64}(1.0e-9)
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

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρE = q[4]
    ρw = q[5]
    Bx = q[6]
    By = q[7]
    Bz = q[8]
    ψ  = q[9]

    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ

    c_h = c_h_mhd[]

    p    = max(pressure_mhd(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz, ψ), p_floor_mhd[])
    magp = 0.5*(Bx*Bx + By*By + Bz*Bz)   # magnetic pressure ½|B|²
    ptot = p + magp                       # total (gas + magnetic) pressure
    vdB  = u*Bx + v*By + w*Bz             # v·B

    # Energy-flux prefactor:
    #   ½ρ|v|² + γp/(γ-1) + |B|²   (= ρE + p + ½|B|² - ½ψ²)
    ke   = 0.5*(ρu*u + ρv*v + ρw*w)
    enfl = ke + γ_mhd*p/(γ_mhd - 1.0) + 2.0*magp

    # x-direction flux
    F[1] = ρu
    F[2] = ρu*u + ptot - Bx*Bx
    F[3] = ρv*u - Bx*By
    F[4] = u*enfl - Bx*vdB + c_h*ψ*Bx
    F[5] = ρw*u - Bx*Bz
    F[6] = c_h*ψ
    F[7] = u*By - v*Bx
    F[8] = u*Bz - w*Bx
    F[9] = c_h*Bx

    # y-direction flux (the paper's z-direction)
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
