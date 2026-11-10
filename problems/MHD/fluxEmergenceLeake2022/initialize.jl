#---------------------------------------------------------------------------------
# Initial state of the 3D flux-emergence parameter study of
#
#   J. E. Leake, M. G. Linton & S. K. Antiochos, "The role of reconnection in the
#   onset of solar eruptions", ApJ 934:10 (2022), Sec. 2.4
#
# in the units of user_flux.jl (L₀ = 170 km, B₀ = 0.13 T, g₀ = 274 m s⁻²;
# p = ρT/m_f with m_f = 1 (below), g = 1, magnetic pressure ½|B|²).
#
# 1. Atmosphere (Eqs. 23-25, as Leake et al. 2013, 2014), hydrostatic, ρ(0) = 1:
#      T = 1 − (γ−1)/γ · z      z ≤ 0     (convection zone, polytropic)
#      T = 1                    0 < z < 10 (photosphere/chromosphere)
#      T = 150^((z−10)/10)      10 ≤ z < 20 (transition region)
#      T = 150                  z ≥ 20    (corona)
#    dp/dz = −ρ = −m_f p/T in closed form per layer (lk22_pressure). m_f (:lk22_mf) is
#    1: the text gives μ_m = 1.25 m_p with T₀ = m_p L₀g₀/k_B, i.e. p = ρT/1.25, but the
#    paper's Fig. 1 (ρ = 1.2e-2 kg m⁻³ at −5.1 Mm, 1e-11 at the base of the corona and
#    a 25 Mm coronal scale height) needs p = ρT, as Leake et al. (2014)'s Fig. 1 does;
#    the reduced mass sits in the temperature unit (README).
#
# 2. Coronal dipole (Eqs. 26-27): a point dipole of strength d at (0, 0, z_d)
#    below the domain, moment along x,
#      B = d [3 (m̂·r̂) r̂ − m̂]/r³,  r = (x, y, z − z_d),
#    d = 15000 (the paper's d = 15,000 B₀L₀³/μ₀) and z_d = −100 (as in Leake et
#    al. 2013, 2014; the paper does not print it). These give |B| = d/z_d³ =
#    0.015 B₀ = 19.5 G at the photosphere above the dipole, the paper's "20 G";
#    with a literal μ₀/4π the same d would give 1.55 G (README). The moment is
#    m̂ = −x̂, so that the dipole field above the rope points along +x (README:
#    the sign that makes θ = 0 the stable, parallel case of the paper's Fig. 2).
#
# 3. Flux rope (Eqs. 28-34), in the frame rotated about z by θ = nπ/4,
#      x_r = cos θ x + sin θ y,   y_r = −sin θ x + cos θ y,
#    along y_r at x_r = 0, z = z_rope = −12, of width w = 2.5 and strength B_t = 5,
#      B_axial = B_t e^(−r²/w²) ŷ_r,  B_twist = (r/w) B_t e^(−r²/w²) (azimuthal),
#    r² = x_r² + (z − z_rope)². The twist is right-handed (the paper's q = −1/w in
#    its (r₁, φ₁, y₁) frame): above the axis it points along +x_r, below along −x_r,
#    so for θ = 0 it is parallel to the dipole field above it and for θ = π
#    antiparallel. The pressure deficit p₁ balances the rope's Lorentz force,
#      p₁ = ½ B_t² e^(−2r²/w²) [(w²/2 − r²)/w² − 1],
#    and the buoyancy is set along the rope by λ = 10 (Eq. 34, with the sign that
#    makes the middle of the rope light, as the text and Leake et al. 2013, 2014):
#      p = p₀(z) + p₁,   ρ = ρ₀(z) [1 + (p₁/p₀(z)) e^(−(y_r/λ)²)].
#
# The reference state qe (well balance, DynSGS departure, positivity floors) is
# the atmosphere plus the dipole, without the rope.
#---------------------------------------------------------------------------------

if !@isdefined(LK22_TCOR)
    const LK22_TCOR = 150.0    # T_cor/T_ph
    const LK22_ZR   = -12.0    # rope axis height z_rope
    const LK22_W    = 2.5      # rope width w
    const LK22_BT   = 5.0      # rope axial field B_t
    const LK22_LAM  = 10.0     # buoyancy length λ
end
if !@isdefined(lk22_d)
    const lk22_d         = Ref{Float64}(15000.0)  # dipole strength (moment along −x)
    const lk22_zd        = Ref{Float64}(-100.0)   # dipole depth z_d
    const lk22_cz_a      = Ref{Float64}(1.0)      # convection-zone gradient factor (1: Eq. 23 as printed)
    const lk22_mf        = Ref{Float64}(1.0)      # m_f of p = ρT/m_f (README)
    const lk22_cθ        = Ref{Float64}(-1.0)     # cos θ, sin θ of the rope frame
    const lk22_sθ        = Ref{Float64}(0.0)
    const lk22_glm_alpha = Ref{Float64}(1.0)      # Dedner damping α_p (rate α_p c_h/Δh)
    const lk22_glm_dh    = Ref{Float64}(1.0)
    const lk22_damp_xd   = Ref{Float64}(207.9)    # side damping from |x|, |y| = 3/8 of the width
    const lk22_damp_zd   = Ref{Float64}(464.2)    # top damping from z = 7/8 max(z)
    const lk22_xmax      = Ref{Float64}(277.2)
    const lk22_zmax      = Ref{Float64}(530.5)
end

# Temperature and pressure of the background atmosphere at height z (item 1).
@inline function lk22_temperature(z)
    a = lk22_cz_a[]
    z <= 0.0  && return 1.0 - a*(γ_mhd - 1.0)/γ_mhd*z
    z < 10.0  && return 1.0
    z < 20.0  && return LK22_TCOR^((z - 10.0)/10.0)
    return LK22_TCOR
end

@inline function lk22_pressure(z)
    a  = lk22_cz_a[]
    p0 = 1.0/lk22_mf[]                                # ρ(0) = 1, T(0) = 1
    if z <= 0.0
        return p0*lk22_temperature(z)^(lk22_mf[]*γ_mhd/(a*(γ_mhd - 1.0)))
    end
    p10 = p0*exp(-lk22_mf[]*10.0)
    z < 10.0 && return p0*exp(-lk22_mf[]*z)
    α = log(LK22_TCOR)/10.0
    z < 20.0 && return p10*exp(-lk22_mf[]*(1.0 - exp(-α*(z - 10.0)))/α)
    p20 = p10*exp(-lk22_mf[]*(1.0 - exp(-α*10.0))/α)
    return p20*exp(-lk22_mf[]*(z - 20.0)/LK22_TCOR)
end

@inline lk22_density(z) = lk22_mf[]*lk22_pressure(z)/lk22_temperature(z)

# Field of the point dipole (item 2): moment −d x̂ at (0, 0, z_d).
@inline function lk22_dipole(x, y, z)
    d  = lk22_d[]
    s  = z - lk22_zd[]
    r2 = x*x + y*y + s*s
    r5 = r2*r2*sqrt(r2)
    return d*(r2 - 3.0*x*x)/r5, -3.0*d*x*y/r5, -3.0*d*x*s/r5
end

# Rope field (B_x, B_y, B_z), pressure deficit p₁ and the along-rope coordinate y_r (item 3).
@inline function lk22_rope(x, y, z)
    c, s = lk22_cθ[], lk22_sθ[]
    xr =  c*x + s*y
    yr = -s*x + c*y
    ζ  = z - LK22_ZR
    r2 = xr*xr + ζ*ζ
    w2 = LK22_W*LK22_W
    Ba = LK22_BT*exp(-r2/w2)
    Bxr = ζ*Ba/LK22_W
    Bz  = -xr*Ba/LK22_W
    p1  = 0.5*LK22_BT^2*exp(-2.0*r2/w2)*((0.5*w2 - r2)/w2 - 1.0)
    return c*Bxr - s*Ba, s*Bxr + c*Ba, Bz, p1, yr
end

# (cos θ, sin θ) of θ = nπ/4, exact for integer n.
function lk22_rotation(n)
    h = sqrt(0.5)
    if isinteger(n)
        return ((1.0, 0.0), (h, h), (0.0, 1.0), (-h, h), (-1.0, 0.0), (-h, -h), (0.0, -1.0), (h, -h))[mod(Int(n), 8) + 1]
    end
    θ = n*π/4
    return cos(θ), sin(θ)
end

function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && println(" Initialize fields for 3D GLM-MHD flux emergence (Leake et al. 2022) ........... ")

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "β", "Bmag", "log10ρ"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    inputs[:backend] != CPU()         && error(" problems/MHD/fluxEmergenceLeake2022: only the CPU backend is supported.")
    inputs[:SOL_VARS_TYPE] != TOTAL() && error(" problems/MHD/fluxEmergenceLeake2022: only SOL_VARS_TYPE = TOTAL() is supported.")

    nsim = get(inputs, :lk22_sim, 4)
    lk22_cθ[], lk22_sθ[] = lk22_rotation(nsim)
    ldip            = get(inputs, :lk22_dipole, true)
    lk22_d[]        = ldip ? Float64(get(inputs, :lk22_d, 15000.0)) : 0.0
    lk22_zd[]       = Float64(get(inputs, :lk22_zd, -100.0))
    lk22_cz_a[]     = Float64(get(inputs, :lk22_cz_gradient, 1.0))
    lk22_mf[]       = Float64(get(inputs, :lk22_mf, 1.0))
    lk22_glm_alpha[] = Float64(get(inputs, :lk22_glm_alpha, 1.0))
    lk22_xmax[]     = MPI.Allreduce(maximum(abs, @view(mesh.coords[1,1:mesh.npoin])), MPI.MAX, comm)
    zmin            = MPI.Allreduce(minimum(@view(mesh.coords[3,1:mesh.npoin])), MPI.MIN, comm)
    lk22_zmax[]     = MPI.Allreduce(maximum(@view(mesh.coords[3,1:mesh.npoin])), MPI.MAX, comm)
    # damping regions (Eqs. 21-22): |x|, |y| > 3L_x/8 (L_x the full width, i.e. the
    # outer eighth of each side) and z > 7 max(z)/8
    lk22_damp_xd[]  = Float64(get(inputs, :lk22_damp_x_frac, 0.75))*lk22_xmax[]
    lk22_damp_zd[]  = Float64(get(inputs, :lk22_damp_z_frac, 7.0/8.0))*lk22_zmax[]
    lk22_zd[] >= zmin && error(" problems/MHD/fluxEmergenceLeake2022: the dipole (z_d = $(lk22_zd[])) must lie below the domain (z ≥ $(zmin)).")
    lrope = get(inputs, :lk22_rope, true)

    γm1 = γ_mhd - 1.0
    ch_local = 0.0
    for ip = 1:mesh.npoin
        x, y, z = mesh.coords[1,ip], mesh.coords[2,ip], mesh.coords[3,ip]
        p0 = lk22_pressure(z)
        ρ0 = lk22_density(z)
        Bxd, Byd, Bzd = lk22_dipole(x, y, z)

        # reference state: atmosphere + dipole
        q.qe[ip,1] = ρ0
        q.qe[ip,2] = 0.0; q.qe[ip,3] = 0.0; q.qe[ip,5] = 0.0
        q.qe[ip,4] = p0/γm1 + 0.5*(Bxd*Bxd + Byd*Byd + Bzd*Bzd)
        q.qe[ip,6] = Bxd; q.qe[ip,7] = Byd; q.qe[ip,8] = Bzd
        q.qe[ip,9] = 0.0
        q.qe[ip,end] = p0

        Bxr, Byr, Bzr, p1, yr = lrope ? lk22_rope(x, y, z) : (0.0, 0.0, 0.0, 0.0, 0.0)
        p  = p0 + p1
        ρ  = ρ0*(1.0 + (p1/p0)*exp(-(yr/LK22_LAM)^2))
        Bx = Bxd + Bxr; By = Byd + Byr; Bz = Bzd + Bzr
        q.qn[ip,1] = ρ
        q.qn[ip,2] = 0.0; q.qn[ip,3] = 0.0; q.qn[ip,5] = 0.0
        q.qn[ip,4] = p/γm1 + 0.5*(Bx*Bx + By*By + Bz*Bz)
        q.qn[ip,6] = Bx; q.qn[ip,7] = By; q.qn[ip,8] = Bz
        q.qn[ip,9] = 0.0
        q.qn[ip,end] = p

        ch_local = max(ch_local, sqrt((γ_mhd*p + Bx*Bx + By*By + Bz*Bz)/ρ))   # fast-speed bound, at rest
    end
    # GLM cleaning speed: the fastest wave of the initial state (the Alfvén speed
    # just above the transition region), held constant.
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    # Smallest nodal spacing (Dedner's damping rate α_p c_h/Δh).
    dh_local = Inf
    for ie = 1:mesh.nelem
        a = mesh.connijk[ie,1,1,1]
        for b in (mesh.connijk[ie,2,1,1], mesh.connijk[ie,1,2,1], mesh.connijk[ie,1,1,2])
            dh_local = min(dh_local, sqrt(sum(abs2, @view(mesh.coords[:,b]) .- @view(mesh.coords[:,a]))))
        end
    end
    lk22_glm_dh[] = MPI.Allreduce(dh_local, MPI.MIN, comm)

    get(inputs, :lwell_balanced, false) || rank == 0 &&
        println(" WARNING problems/MHD/fluxEmergenceLeake2022: :lwell_balanced is off; the discrete",
                " imbalance of the atmosphere and the dipole will drive flows from t = 0.")

    if rank == 0
        Bph = ldip ? lk22_d[]/abs(lk22_zd[])^3 : 0.0
        println(" Simulation $(nsim): θ = $(nsim)π/4 (cos θ = $(lk22_cθ[]), sin θ = $(lk22_sθ[])); rope $(lrope ? "on" : "OFF (:lk22_rope => false)")")
        println(" Dipole: $(ldip ? "d = $(lk22_d[]) (moment along −x) at z_d = $(lk22_zd[]); |B| at (0,0,0) = $(Bph) B₀ = $(round(1300*Bph, digits=2)) G" : "OFF (:lk22_dipole => false)")")
        println(" Atmosphere (m_f = $(lk22_mf[])): p(-30) = $(lk22_pressure(-30.0)), ρ(-30) = $(lk22_density(-30.0)); ρ(20) = $(lk22_density(20.0)), ρ(top) = $(lk22_density(lk22_zmax[]))")
        println(" Rope: B_t = $(LK22_BT), w = $(LK22_W), z_rope = $(LK22_ZR), λ = $(LK22_LAM)")
        println(" GLM c_h = $(c_h_mhd[]); smallest nodal spacing Δh = $(lk22_glm_dh[]); ψ damping rate α_p c_h/Δh = $(lk22_glm_alpha[]*c_h_mhd[]/lk22_glm_dh[]) per t₀")
        println(" Velocity damping: |x|, |y| > $(lk22_damp_xd[]) (boundary at $(lk22_xmax[])) and z > $(lk22_damp_zd[]) (top at $(lk22_zmax[]))")
        println(" Initialize fields for 3D GLM-MHD flux emergence (Leake et al. 2022) ........... DONE")
        flush(stdout)
    end
    return q
end
