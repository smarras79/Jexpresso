#---------------------------------------------------------------------------------
# 2D ideal GLM-MHD equations
#
# The ideal MHD system solved by
#
#   A. Bormanis, C. A. Leon, A. Scheinker,
#   "Solving the Orszag-Tang vortex magnetohydrodynamics problem with
#    physics-constrained convolutional neural networks",
#   Phys. Plasmas 31, 012101 (2024), Eqs. (1)-(4) and Appendix B, Eqs. (B3)-(B6),
#
# augmented with Dedner's generalized-Lagrange-multiplier (GLM) divergence
# cleaning. The paper keeps ∇·B = 0 with constrained transport on a
# staggered finite-volume mesh (its Sec. III A); that machinery has no
# analogue in a collocated continuous-Galerkin spectral element method, so
# here the divergence constraint is enforced the way the DG/SEM MHD
# literature does it — by an extra transported field ψ that radiates
# divergence error out of the domain at the constant speed c_h. This is
# the same equation set already used by
# problems/MHD/kelvinHelmholtzChan2022.
#
# State (2D in space, but the third components of velocity and magnetic
# field are kept because plasmas admit 3D electromagnetic interactions in
# 2D problems):
#
#   q = (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)
#
# NOTICE on the ordering: the MHD literature writes the state as
# (ρ, ρv, E, B, ψ) with the three momentum components contiguous. Here the
# total energy ρE sits in slot 4 and the out-of-plane momentum ρw in slot 5
# because Jexpresso's shared 2D kernels (viscous τ·u augmentation and the
# sound-speed/CFL diagnostic in soundSpeed.jl) assume that slot 4 of a 2D
# system carries the energy. Slots 5..9 are treated by the SGS machinery as
# generic advected scalars, which is exactly what we want for ρw, B and ψ
# in this simple Smagorinsky setup.
#
# For the Orszag-Tang vortex w ≡ Bz ≡ 0 for all time (nothing in the 2D
# fluxes couples them once they vanish), so those two slots simply stay at
# zero; they are carried for generality of the equation set.
#
# Pressure:  p = (γ-1)·(ρE - ½ρ|v|² - ½|B|² - ½ψ²)
#
# This is the CONSERVATIVE part of the GLM-MHD system only. The
# non-conservative Powell/Galilean-GLM term is required for entropy
# stability of split-form/two-point-flux discretizations, which are NOT
# used here (no ES/KEP; plain weak form + Smagorinsky SGS dissipation). It
# is therefore omitted, exactly as in the companion KHI case.
#---------------------------------------------------------------------------------

# Heat-capacity ratio. Bormanis et al. take a monatomic perfect gas with
# γ = 5/3 (their Sec. I A, "we ... let γ = 5/3 to have better analogy to
# the 3D case"). The Orszag-Tang initial condition in initialize.jl is
# written in terms of γ, so changing this value changes ρ and p too.
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
