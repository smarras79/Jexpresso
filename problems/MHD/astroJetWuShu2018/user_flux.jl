#---------------------------------------------------------------------------------
# 2D ideal GLM-MHD equations — magnetized astrophysical jet.
#
# Same equation set as problems/MHD/orszagTangBormanis2024 and
# problems/MHD/fluxEmergenceSon2025 (ideal MHD + Dedner's generalized-Lagrange-
# multiplier divergence cleaning), specialized to γ = 1.4. See
# problems/MHD/orszagTangBormanis2024/EQUATIONS.md for the derivation.
#
# State (2D in space; the third components of v and B are carried because a
# 2D plasma still admits out-of-plane electromagnetic interaction):
#
#   q = (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)
#
# NOTICE on the ordering: the MHD literature writes (ρ, ρv, E, B, ψ) with the
# three momentum components contiguous. Here the total energy ρE sits in slot
# 4 and the out-of-plane momentum ρw in slot 5, because Jexpresso's shared 2D
# kernels (the viscous τ·u augmentation in rhs.jl and the sound-speed/CFL
# diagnostic in soundSpeed.jl) assume slot 4 of a 2D system carries the energy.
#
# Pressure:  p = (γ-1)·(ρE - ½ρ|v|² - ½|B|² - ½ψ²)
#
# This is the CONSERVATIVE part of the GLM-MHD system only. The
# non-conservative Powell/Galilean-GLM term is needed for entropy stability of
# split-form/two-point-flux discretizations, which are NOT used here (no
# ES/KEP; plain weak form + DynSGS dissipation), so it is omitted exactly as
# in the three companion MHD cases.
#---------------------------------------------------------------------------------

# Heat-capacity ratio. Wu & Shu (SISC 40(5):B1302, 2018, Example 5.6) set
# γ = 1.4 for this test, and the whole initial condition is written in terms of
# it: the ambient density is 0.1γ and the beam density is γ, which is what
# makes the beam sound speed sqrt(γ p / ρ) = sqrt(γ·1/γ) EXACTLY 1 and the
# injection speed 800 exactly Mach 800. Changing this constant changes the
# problem, not just the gas.
if !@isdefined(γ_mhd)
    const γ_mhd = 1.4
end
if abs(γ_mhd - 1.4) > 1e-12
    error(" problems/MHD/astroJetWuShu2018: γ_mhd = $(γ_mhd) is already defined by another MHD case loaded in this Julia session (this case needs γ = 1.4). Restart Julia before running this case.")
end

# Hyperbolic divergence-cleaning speed c_h. Filled by initialize.jl with the
# maximum wave speed |v| + c_f over the initial condition AND the prescribed
# nozzle inflow state, then kept constant for the whole run (standard GLM
# practice). The inflow has to be included: the domain is at rest at t = 0, so
# the initial condition alone would give c_h ≈ 38 (B_a = √200) while the beam
# injected through the nozzle from the first step travels at 800.
if !@isdefined(c_h_mhd)
    const c_h_mhd = Ref{Float64}(1.0)
end

# Pressure floor used ONLY inside the flux evaluation — the same device as
# problems/MHD/fluxEmergenceSon2025, and this case needs it more than any
# other in the tree.
#
# In the beam the gas pressure is a 5.6 PARTS PER MILLION residue of the total
# energy: with ρ = 1.4, |v| = 800, p = 1, B_a = √200,
#
#   ρE = p/(γ-1) + ½ρ|v|² + ½|B|² = 2.5 + 448000 + 100 = 448102.5
#
# so p = 0.4·(ρE − 448100), and a RELATIVE error of 5.6·10⁻⁶ in ρE wipes the
# pressure out entirely — as does 2.8·10⁻⁶ in the momentum ρv, since the kinetic
# energy responds to it with a factor 2. That is the whole difficulty of
# this benchmark and the reason the papers that run it use an explicit
# positivity-preserving limiter. There is none here for MHD: the node-wise
# realizability repair of src/kernel/positivity/ (:lpositivity) is scoped to
# neqs == nsd + 2 exactly and errors on this nine-field state, because it knows
# nothing about the magnetic energy that p = (γ-1)(ρE − KE − ½|B|² − ½ψ²)
# subtracts. DynSGS is therefore the only stabilization, and this floor exists
# so that a momentarily inadmissible node still produces a FINITE flux instead
# of poisoning the whole RHS. It is not a
# positivity guarantee and it is not conservative where it fires — if it fires
# at all, the run is already in trouble and the output must not be trusted;
# JEXPRESSO_AJ_PFLOOR=0 disables it to see the unshielded behaviour.
#
# 10⁻⁶ is 10⁻⁶ of the (uniform) initial pressure, far below anything the
# problem produces physically: the paper's log₁₀p figures bottom out near
# p ≈ 10⁻¹ in the rarefied cocoon.
if !@isdefined(p_floor_mhd)
    const p_floor_mhd = Ref{Float64}(
        something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_PFLOOR", "")), 1.0e-6))
end

# Density guard for the DIVISIONS in the flux, u = ρu/ρ and ke = |ρv|²/(2ρ).
#
# With :lpositivity => true (this case's default) this is unreachable: the
# realizability repair runs in rhs! before any flux is evaluated and leaves
# ρ ≥ :positivity_rho_min = 1.4e-7 everywhere, six orders above this. It exists
# for a run with the repair turned off, and it exists because an UNGUARDED
# division is how a single bad node stops being a local defect: ρ → 0 gives
# u = ±Inf, the flux goes Inf, and with :dsgs_norms => "domain" the very next
# DynSGS reduction (an Allreduce over the whole domain) turns that Inf into a NaN
# in ⟨q⟩ and hence in ν on EVERY element of EVERY rank. That is the mechanism
# behind "100.0% of local entries non-finite" on all ranks in one step: the
# failure is local, the reduction makes it global. This floor keeps the flux
# finite so the failure stays local and locatable. It is not a positivity
# guarantee and it is not conservative where it fires — the repair is.
if !@isdefined(ρ_floor_mhd)
    const ρ_floor_mhd = Ref{Float64}(
        something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_RHOFLOOR", "")), 1.0e-14))
end

@inline function pressure_mhd(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz, ψ)
    γm1 = γ_mhd - 1.0
    ke  = 0.5*(ρu*ρu + ρv*ρv + ρw*ρw)/max(ρ, ρ_floor_mhd[])
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

    # guarded: see ρ_floor_mhd above. Unreachable with :lpositivity => true.
    ρg = max(ρ, ρ_floor_mhd[])
    u = ρu/ρg
    v = ρv/ρg
    w = ρw/ρg

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

#---------------------------------------------------------------------------------
# THE PROBLEM'S CONSTANTS, in one place.
#
# Wu & Shu, SIAM J. Sci. Comput. 40(5):B1302-B1329 (2018), Example 5.6
# ("Astrophysical jets"), which adds a magnetic field to the Mach 800 gas
# dynamical jet of Balsara, JCP 231:7504-7517 (2012):
#
#   "Initially, the physical domain [-0.5,0.5]x[0,1.5] is filled with a uniform
#    static medium with density 0.1γ and unit pressure, and the adiabatic index
#    γ is set as 1.4. Through the inlet part (|x| < 0.05) on the bottom
#    boundary (y = 0), a dense jet with speed 800 is injected in the
#    y-direction with a density of γ and a pressure equal to the ambient
#    pressure. The fixed inflow beam condition is specified on the nozzle
#    {y = 0, |x| < 0.05}, and the others are outflow boundary conditions. We
#    initialize the magnetic field with magnitude B_a along the y-direction."
#
# Their three configurations, all at the same γ, domain and beam:
#
#   (i)   moderately magnetized        B_a = √200     β_a = 2p/B_a² = 10⁻²
#   (ii)  strongly magnetized          B_a = √2000    β_a = 10⁻³
#   (iii) extremely strongly magnetized B_a = √20000  β_a = 10⁻⁴
#
# Note how the numbers are built out of γ. The beam sound speed is
# sqrt(γ p/ρ) = sqrt(γ·1/γ) = 1 EXACTLY, so "speed 800" is exactly Mach 800;
# the ambient sound speed is sqrt(1.4/0.14) = sqrt(10) ≈ 3.1623, so the beam
# is also a factor 10 denser than what it ploughs into.
#
# Overrides (read once, at include time; also settable from the REPL through
# the Refs):
#   JEXPRESSO_AJ_BA2    B_a SQUARED, as the paper quotes it (default 200)
#   JEXPRESSO_AJ_UJET   injection speed = beam Mach number (default 800)
#   JEXPRESSO_AJ_SMOOTH nozzle-lip transition HALF-WIDTH s (default: one
#                       element, resolved from the mesh; 0 = the paper's exact
#                       top hat, which does not run — see user_bc.jl)
#   JEXPRESSO_AJ_TRAMP  inflow turn-on time τ (default 2h/u_jet = 125 steps;
#                       0 = the impulsive start, which does not run either)
#---------------------------------------------------------------------------------
if !@isdefined(AJ_XNOZZLE)
    const AJ_XNOZZLE = 0.05             # nozzle half-width
    const AJ_RHO_AMB = 0.1*γ_mhd        # = 0.14
    const AJ_P_AMB   = 1.0              # "unit pressure"
    const AJ_RHO_JET = γ_mhd            # = 1.4  → beam sound speed exactly 1
    const AJ_P_JET   = 1.0              # "a pressure equal to the ambient pressure"
end
if !@isdefined(aj_Ba)
    const aj_Ba     = Ref{Float64}(sqrt(something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_BA2",   "")), 200.0)))
    const aj_ujet   = Ref{Float64}(     something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_UJET",  "")), 800.0))
    # -1 is the AUTO sentinel: initialize.jl resolves it to ONE element from the
    # mesh it was actually given, so the transition spans two elements. An explicit
    # JEXPRESSO_AJ_SMOOTH=0 means the paper's exact top hat, which is measured not
    # to run — see user_bc.jl and README.md §10-11.
    const aj_smooth = Ref{Float64}(     something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_SMOOTH", "")), -1.0))
    # Inflow turn-on time. -1 is the AUTO sentinel, resolved by initialize.jl to
    # 2h/u_jet — the time the beam needs to cross two elements, which is 125 time
    # steps on either shipped mesh because Δt scales with h. 0 is the impulsive
    # start of the paper, which is measured not to run: see user_bc.jl.
    const aj_tramp  = Ref{Float64}(     something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_TRAMP",  "")), -1.0))
end

# The two states of the problem, as CONSERVED 9-tuples in this case's slot
# order (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ). Bx = Bz = ψ = 0 and By = B_a in
# both, so ∇·B = ∂_y B_a = 0 identically — the initial field and the injected
# field are divergence free by construction, not to discretization accuracy.
@inline function aj_ambient_state()
    Ba = aj_Ba[]
    ρ  = AJ_RHO_AMB
    ρE = AJ_P_AMB/(γ_mhd - 1.0) + 0.5*Ba*Ba
    return (ρ, 0.0, 0.0, ρE, 0.0, 0.0, Ba, 0.0, 0.0)
end

@inline function aj_jet_state()
    Ba = aj_Ba[]
    uj = aj_ujet[]
    ρ  = AJ_RHO_JET
    ρE = AJ_P_JET/(γ_mhd - 1.0) + 0.5*ρ*uj*uj + 0.5*Ba*Ba
    return (ρ, 0.0, ρ*uj, ρE, 0.0, 0.0, Ba, 0.0, 0.0)
end

# Maximum wave speed |v| + c_f of a conserved state, with the fast
# magnetosonic speed bounded over all propagation directions,
# c_f ≤ sqrt(a² + b²), a² = γp/ρ (sound), b² = |B|²/ρ (Alfvén). This is the
# same bound compute_dsgs_viscosity!(::DSGS_MHD) uses for its wave-speed cap.
@inline function aj_wave_speed(s)
    ρ  = s[1]
    u  = s[2]/ρ; v = s[3]/ρ; w = s[5]/ρ
    B2 = s[6]*s[6] + s[7]*s[7] + s[8]*s[8]
    p  = pressure_mhd(ρ, s[2], s[3], s[5], s[4], s[6], s[7], s[8], s[9])
    return sqrt(u*u + v*v + w*w) + sqrt(max(γ_mhd*p/ρ + B2/ρ, 0.0))
end
