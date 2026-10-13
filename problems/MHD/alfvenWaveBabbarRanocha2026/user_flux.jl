#---------------------------------------------------------------------------------
# 2D ideal GLM-MHD equations: the conservative flux of
# problems/MHD/orszagTangBormanis2024 and problems/MHD/smoothVortex, unchanged.
#
#   q = (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)
#
# The total energy ρE sits in slot 4 and the out-of-plane momentum ρw in slot 5
# because Jexpresso's shared 2D kernels (viscous τ·u augmentation and the
# sound-speed/CFL diagnostic in soundSpeed.jl) assume that slot 4 of a 2D
# system carries the energy.
#
# Unlike the Orszag-Tang vortex, the Alfvén wave lives in the out-of-plane
# components as much as in the plane: w = B_z = 0.1 cos(2π(x⊥ + t)), so the
# F[5], F[8], G[5], G[8] entries below carry half of the wave.
#
# Pressure:  p = (γ-1)·(ρE - ½ρ|v|² - ½|B|² - ½ψ²)
#
# Babbar & Ranocha solve the GLM-MHD system with its non-conservative part
# (their eq. 81: the Godunov-Powell term ∝ ∇·B and the GLM transport v·∇ψ).
# Both terms vanish on this exact solution (∇·B = 0, ψ = 0), so the
# conservative system below has the same exact solution; they differ only in
# how the discretization's own ∇·B and ψ errors are carried, which is what the
# ψ damping of user_source.jl takes care of here.
#---------------------------------------------------------------------------------

# Heat-capacity ratio, γ = 5/3 for this test.
if !@isdefined(γ_mhd)
    const γ_mhd = 5.0/3.0
end

# Hyperbolic divergence-cleaning speed c_h. Filled by initialize.jl with the
# maximum wave speed of the initial condition and kept constant throughout
# the run (standard GLM practice).
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
