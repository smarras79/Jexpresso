#---------------------------------------------------------------------------------
# 2D ideal GLM-MHD equations
#
# Chan, Ranocha, Rueda-Ramírez, Gassner, Warburton,
# "On the entropy projection and the robustness of high order entropy stable
#  discontinuous Galerkin schemes for under-resolved flows",
# Frontiers in Physics 10:898028 (2022), Section 3.2, Eq. (6).
#
# State (2D in space, but the third components of velocity and magnetic field
# are kept because plasmas admit 3D electromagnetic interactions in 2D
# problems):
#
#   q = (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)
#
# NOTICE on the ordering: the paper writes the state as
# (ρ, ρv, E, B, ψ) with the three momentum components contiguous. Here the
# total energy ρE sits in slot 4 and the out-of-plane momentum ρw in slot 5
# because Jexpresso's shared 2D kernels (viscous τ·u augmentation and the
# sound-speed/CFL diagnostic in soundSpeed.jl) assume that slot 4 of a 2D
# system carries the energy. Slots 5..9 are treated by the SGS machinery as
# generic advected scalars, which is exactly what we want for ρw, B and ψ
# in this simple Smagorinsky setup.
#
# ψ is the generalized Lagrange multiplier (GLM) divergence-cleaning field
# and c_h the (constant) hyperbolic divergence-cleaning speed, chosen as the
# maximum wave speed of the initial condition (set once in initialize.jl,
# see the paper's Section 3.2.1).
#
# Pressure:  p = (γ-1)·(ρE - ½ρ|v|² - ½|B|² - ½ψ²)
#
# This is the CONSERVATIVE part of the GLM-MHD system only. The
# non-conservative term Υ of the paper's Eq. (7) — the Powell term
# (∇·B)(0, B, v·B, v, 0) plus the Galilean GLM term (0, 0, ψ v·∇ψ, 0, v·∇ψ) —
# is required for entropy stability of the split-form/two-point-flux
# discretizations, which are NOT used here (no ES/KEP; plain weak form +
# Smagorinsky SGS dissipation). It is therefore omitted for now.
#---------------------------------------------------------------------------------

# Heat-capacity ratio. The paper does not state a numerical value for γ;
# 5/3 is the standard ideal (monatomic) MHD choice used by the GLM-MHD
# literature the paper builds on. Set to 1.4 to match the companion
# CompEuler/kelvinHelmholtzChan2022 Euler case.
if !@isdefined(γ_mhd)
    const γ_mhd = 5.0/3.0
end

# Hyperbolic divergence-cleaning speed c_h. Filled by initialize.jl with the
# maximum wave speed of the initial condition and kept constant throughout
# the run (standard GLM practice, see Section 3.2.1 of the paper).
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

    p    = pressure_mhd(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz, ψ)
    magp = 0.5*(Bx*Bx + By*By + Bz*Bz)   # magnetic pressure ½|B|²
    ptot = p + magp                       # total (gas + magnetic) pressure
    vdB  = u*Bx + v*By + w*Bz             # v·B

    # Energy-flux prefactor of the paper's Eq. (6):
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

    # y-direction flux
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
