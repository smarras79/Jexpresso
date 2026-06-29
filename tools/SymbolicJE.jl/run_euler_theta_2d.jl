#=---------------------------------------------------------------------------------
# run_euler_theta_2d.jl
#
# Stand-alone 2D example for SymbolicJE.jl — the *system* analogue of
# run_advdiff_2d.jl. It writes the 2D compressible Euler equations in
# potential-temperature (θ) form as four coupled conservation laws in unicode and
# time-marches the classic rising thermal bubble of problems/CompEuler/theta.
#
#   ∂ρ /∂t + ∂x(ρu)            + ∂y(ρv)            = 0
#   ∂ρu/∂t + ∂x(ρu·u + p′)     + ∂y(ρu·v)          = 0
#   ∂ρv/∂t + ∂x(ρv·u)          + ∂y(ρv·v + p′)     = -ρ′ g
#   ∂ρθ/∂t + ∂x(ρθ·u)          + ∂y(ρθ·v)          = 0
#
#       u = ρu/(ρ̄+ρ′),  v = ρv/(ρ̄+ρ′),  p = C₀ (ρ̄θ̄+ρθ)^γ,  p′ = p − p̄
#
# As in Jexpresso's `theta` (`:SOL_VARS_TYPE => PERT()`) we evolve the PERTURBATION
# about the hydrostatic background ρ̄(y), (ρθ)‾(y), p̄(y) (∂p̄/∂y = −ρ̄g). Carrying
# the pressure perturbation p′ and the buoyancy source −ρ′g makes the scheme
# DISCRETELY well balanced: the resting background is preserved to machine
# precision (see runtests.jl). A small constant viscosity ν regularizes the
# momentum/θ equations, exactly the role of Jexpresso's `:lvisc => true, :μ`.
#
# Run from the Jexpresso root (so Plots resolves through the project):
#       julia --project=. tools/SymbolicJE.jl/run_euler_theta_2d.jl
#
# Defaults to the self-contained structured finite-difference backend (no gmsh /
# no Jexpresso needed). To run instead on a spectral-element gmsh grid, set
# :method => :sem and :gmsh_filename (see the commented block at the bottom).
#---------------------------------------------------------------------------------=#

isdefined(Main, :SymbolicJE) || include(joinpath(@__DIR__, "src", "SymbolicJE.jl"))
using .SymbolicJE

#---------------------------------------------------------------------------------
# 1. Thermodynamics (PhysicalConst defaults; see globalConstantsPhysics.jl).
#    γ, C₀, g are ordinary numbers (the exponent in (·)^γ and the buoyancy factor);
#    only the field-valued quantities are symbols.
#---------------------------------------------------------------------------------
const Rair = 287.0; const cp = 1004.0; const cv = 717.0
const γ    = cp / cv
const g    = 9.80616
const pref = 101200.0
const C0   = (Rair^γ) / pref^(γ - 1.0)
const θref = 300.0

# hydrostatic background with constant θ = θref  (∂p̄/∂y = −ρ̄g satisfied analytically)
p̄(x, y)  = pref * (1.0 - g * y / (cp * θref))^(cp / Rair)
ρθ̄(x, y) = (p̄(x, y) / C0)^(1 / γ)          # = ρ̄ θref
ρ̄(x, y)  = ρθ̄(x, y) / θref

#---------------------------------------------------------------------------------
# 2. The equations, with live Julia symbols (residual form, = 0 implied). The
#    unknowns are the PERTURBATIONS ρ, ρu, ρv, ρθ; ρb, ρθb, pb are the background
#    fields (resolved from inputs as functions of (x, y)); ν is the viscosity.
#---------------------------------------------------------------------------------
@vars ρ ρu ρv ρθ ρb ρθb pb ν

u  = ρu / (ρb + ρ)                         # velocity from total momentum / total density
v  = ρv / (ρb + ρ)
p′ = C0 * (ρθb + ρθ)^γ - pb                # pressure perturbation p − p̄

eqρ  = ∂t(ρ)  + ∂x(ρu)             + ∂y(ρv)
eqρu = ∂t(ρu) + ∂x(ρu*u + p′)      + ∂y(ρu*v)          - ν*Δ(ρu)
eqρv = ∂t(ρv) + ∂x(ρv*u)           + ∂y(ρv*v + p′)     - ν*Δ(ρv)  + ρ*g
eqρθ = ∂t(ρθ) + ∂x((ρθb + ρθ)*u)   + ∂y((ρθb + ρθ)*v)  - ν*Δ(ρθ)

equations = [eqρ, eqρu, eqρv, eqρθ]

#---------------------------------------------------------------------------------
# 3. Inputs — the rising thermal bubble (problems/CompEuler/theta initialize.jl):
#    warm anomaly Δθ = θc(1 − r/r0) inside r0 of (xc, yc), at rest, on a
#    hydrostatic background. We store the PERTURBATIONS (ρ′, ρu, ρv, (ρθ)′).
#---------------------------------------------------------------------------------
function user_inputs()
    xc, yc, r0, θc = 5000.0, 2500.0, 2000.0, 2.0      # bubble (theta/initialize.jl)

    function q0(x, y)
        r  = sqrt((x - xc)^2 + (y - yc)^2)
        Δθ = r < r0 ? θc * (1.0 - r / r0) : 0.0
        θ  = θref + Δθ
        p  = pref * (1.0 - g * y / (cp * θ))^(cp / Rair)
        ρ  = (p / C0)^(1 / γ) / θ                     # total density from (θ, p)
        ρr = ρ̄(x, y)                                  # background density at rest
        return (ρ - ρr, 0.0, 0.0, ρ * θ - ρr * θref)  # (ρ′, ρu, ρv, (ρθ)′)
    end

    # diagnostic: potential-temperature perturbation Δθ = (ρ̄θ̄+ρθ)/(ρ̄+ρ′) − θref
    diag(Q, m) = [ (ρθ̄(m.x[ip], m.y[ip]) + Q[4][ip]) /
                   (ρ̄(m.x[ip], m.y[ip]) + Q[1][ip]) - θref  for ip in 1:m.npoin ]

    return Dict(
        #-------------------------------------------------------------------
        # Space, discretization, structured grid over the bubble box
        #-------------------------------------------------------------------
        :nsd        => 2,
        :method     => :sem,                    # structured FD (self-contained)
        :periodic   => false,
        :xmin       => 0.0,  :xmax => 10000.0,
        :ymin       => 0.0,  :ymax => 10000.0,
        :npoinx     => 100,  :npoiny => 100,
        #-------------------------------------------------------------------
        # Background fields + viscosity referenced by the equations
        #-------------------------------------------------------------------
        :ρb         => ρ̄,
        :ρθb        => ρθ̄,
        :pb         => p̄,
        :ν          => 15.0,                   # constant artificial viscosity (m²/s)
        #-------------------------------------------------------------------
        # Initial condition and the plotted/printed diagnostic
        #-------------------------------------------------------------------
        :q0         => q0,
        :diag       => diag,
        :diag_name  => "dtheta",
        #-------------------------------------------------------------------
        # Time integration (acoustic CFL ⇒ Δt; raise :tend to ~700 s for the
        # full bubble rise — kept short here so the example finishes quickly)
        #-------------------------------------------------------------------
        :tend       => 1000.0,
        :wave_speed => 360.0,                  # ≈ sound speed for the CFL Δt
        :CFL        => 0.5,
        #-------------------------------------------------------------------
        # Output (PNG node-map snapshots of Δθ + on-the-fly window)
        #-------------------------------------------------------------------
        :output_dir => joinpath(@__DIR__, "output_euler_theta"),
        :outformat  => "png",
        :ndiagnostics_outputs => 100,
        :plot_live  => true,
    )
end

# --- spectral-element alternatives (Jexpresso's SEM infrastructure) -------------
# The directional operators ∂x/∂y (and Δ) dispatch on the mesh's discretization,
# so switching the backend is purely a matter of the grid keys — the equations
# above are unchanged. On a SEM mesh, ∂x/∂y become Jexpresso's weak-form
# derivative (rhs.jl's _expansion_inviscid! + DSS + Minv), exactly as in 1D, and
# ∂x(Fx)+∂y(Fy) reproduces the fused inviscid divergence by linearity.
#
#   (a) structured tensor-product SEM (reuses the SAME 1D Jexpresso LGL basis,
#       no gmsh file needed) — swap the :npoinx/:npoiny keys above for:
#         :method => :sem, :nop => 4, :nelx => 25, :nely => 25,
#
#   (b) a gmsh grid read through Jexpresso's sem_setup (full curvilinear metric
#       terms dξdx,…,Je on the EXISTING mesh) — add:
#         :method => :sem, :nop => 4, :lread_gmsh => true,
#         :gmsh_filename => "./meshes/gmsh_grids/hexa_TFI_10x10.msh",
#       (this is the grid problems/CompEuler/theta itself runs on).

#---------------------------------------------------------------------------------
# 4. Solve the system.
#---------------------------------------------------------------------------------
mesh, Q0, Q = SymbolicJE.solve(equations, user_inputs())
