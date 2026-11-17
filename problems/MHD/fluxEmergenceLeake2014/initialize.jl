#---------------------------------------------------------------------------------
# Initial state of the 3D flux-emergence simulations of
#
#   J. E. Leake, M. G. Linton & S. K. Antiochos, "Simulations of emerging magnetic
#   flux. II. The formation of unstable coronal flux ropes and the initiation of
#   coronal mass ejections", ApJ 787:46 (2014)                        (Paper II)
#
# whose equations, domain, boundaries and initial state are "exactly the same as
# used in Paper I" (J. E. Leake, M. G. Linton & T. Török, ApJ 778:99, 2013,
# Sec. 2.4), except that the dipole is reversed. Units of user_flux.jl (L₀, B₀,
# g₀; p = ρT/m_f with m_f = 1 (below), g = 1, magnetic pressure ½|B|²).
#
# 1. Atmosphere (Paper I Eqs. 11-14, hydrostatic, ρ(0) = 1):
#      T = 1 − a (γ−1)/γ · z        z ≤ 0     (convection zone, a = 1)
#      T = 1                        0 < z < 10 (photosphere/chromosphere)
#      T = 150^((z−10)/10)          10 ≤ z < 20 (transition region)
#      T = 150                      z ≥ 20    (corona)
#    dp/dz = −ρ = −m_f p/T integrates in closed form in every layer:
#      z ≤ 0:  p = p₀ T^(m_f γ/(a(γ−1)))   (a polytrope), p₀ = 1/m_f
#      then exponentials (lk_atmosphere below).
#    m_f (:lk_mf) is 1: the papers' text gives μ_m = 1.25 m_p with T₀ = m_p L₀g₀/k_B,
#    i.e. p = ρT/1.25, but their Fig. 1 (pressure 609 at z = −30 and 6e-6 in the
#    corona) and Paper I's statement that a = 1 is exactly adiabatic both need
#    p = ρT: the reduced mass sits in their temperature unit (README).
#    :lk_mf => 1.25 gives the literal-text atmosphere (20 times less dense corona).
#
# 2. Dipole (Paper I Eq. 15, invariant along y): B = ∇×(A_y e_y),
#      A_y = B_d (z − z_d)/r₁³,  r₁² = x² + (z − z_d)²,  z_d = −100,
#    i.e. B_x = B_d (3s² − r₁²)/r₁⁵, B_z = −3 B_d s x/r₁⁵ with s = z − z_d.
#    Paper II reverses it (B_x < 0 above the tube, antiparallel to the tube's top
#    field) and sets its strength by the azimuthal-flux ratio Φ_dip/Φ_tube =
#    0.13, 0.1, 0.067 (SD, MD, WD; its Eqs. 9-10). B_d is solved for here from that
#    ratio on the x = y = 0 line (lk_solve_dipole), with Φ_tube the azimuthal flux of
#    the whole tube, ∫|B_x,tube| dz = B_t R (:lk_flux_def => "whole"): the B_d this
#    gives (16006, 12351, 8314) reproduce the β profiles of Paper II's Fig. 1(a) to
#    3 % at every height, while the printed top-half integral (:lk_flux_def => "top")
#    gives half of them, with β four times the figure's (README).
#
# 3. Flux tube (Paper I Eqs. 16-18), along y at x = 0, z = z_t = −12:
#      B_y = B_t e^(−r²/R²),  B_θ = q r B_y,  r² = x² + (z − z_t)²,
#    B_t = 5, R = 2.5, q = 1/R (right-handed: B_θ along e_y × e_r), so
#      B_x = q (z − z_t) B_y,  B_z = −q x B_y  (B_x > 0 at the top of the tube).
#    Pressure balance ∇p₁ = j × B of the tube field:
#      p₁ = ½ B_t² e^(−2r²/R²) [q²(R²/2 − r²) − 1].
#    Buoyancy (Paper I Eqs. 19-20): p = p₀(z) + p₁,
#      ρ = ρ₀(z) [1 + (p₁/p₀(z)) e^(−(y/λ)²)],  λ = 10,
#    so the tube is buoyant at y = 0 (isothermal with its surroundings) and
#    neutrally buoyant at the y boundaries (same density, lower temperature).
#
# The reference state qe (DynSGS departure, positivity floors, sponge) is the
# atmosphere plus the dipole, without the tube.
#
# 4. Discrete equilibrium: the deck sets :lwell_balanced, so the kernel subtracts
#    the element RHS of qe at every stage (src/kernel/operators/rhs.jl,
#    _well_balanced_reference!). Neither the stratification nor the dipole is an
#    exact equilibrium of the discrete operator; in the corona (β ~ 10⁻², ρ ~ 10⁻⁹)
#    the interpolation error of either is a large acceleration.
#---------------------------------------------------------------------------------

if !@isdefined(LK_TCOR)
    const LK_TCOR  = 150.0    # T_cor/T_ph
    const LK_ZT    = -12.0    # tube axis height
    const LK_RT    = 2.5      # tube radius R (Paper II: a)
    const LK_BT    = 5.0      # tube axial field B_t
    const LK_LAM   = 10.0     # buoyancy length λ along the tube
    const LK_ZD    = -100.0   # dipole source height z_d
    const LK_RATIO = Dict("SD" => 0.13, "MD" => 0.1, "WD" => 0.067, "ND" => 0.0)
end
if !@isdefined(lk_Bd)
    const lk_Bd          = Ref{Float64}(0.0)     # signed B_d (set by initialize)
    const lk_cz_a        = Ref{Float64}(1.0)     # a of Paper I Eq. 11
    const lk_mf          = Ref{Float64}(1.0)     # m_f of p = ρT/m_f (README)
    const lk_glm_alpha   = Ref{Float64}(1.0)     # Dedner damping α_p (rate α_p c_h/Δh)
    const lk_glm_dh      = Ref{Float64}(1.0)
    const lk_damp_xd     = Ref{Float64}(96.0)    # side damping starts at |x|, |y| = 96
    const lk_damp_zd     = Ref{Float64}(180.0)   # top damping starts at z = 180 (Paper II; Paper I 170)
    const lk_xmax        = Ref{Float64}(127.6)
    const lk_zmax        = Ref{Float64}(209.8)
end

# Temperature and pressure of the background atmosphere at height z (Section 1 above).
@inline function lk_temperature(z)
    a = lk_cz_a[]
    z <= 0.0  && return 1.0 - a*(γ_mhd - 1.0)/γ_mhd*z
    z < 10.0  && return 1.0
    z < 20.0  && return LK_TCOR^((z - 10.0)/10.0)
    return LK_TCOR
end

@inline function lk_pressure(z)
    a  = lk_cz_a[]
    p0 = 1.0/lk_mf[]                                  # ρ(0) = 1, T(0) = 1
    if z <= 0.0
        return p0*lk_temperature(z)^(lk_mf[]*γ_mhd/(a*(γ_mhd - 1.0)))
    end
    p10 = p0*exp(-lk_mf[]*10.0)
    z < 10.0 && return p0*exp(-lk_mf[]*z)
    α = log(LK_TCOR)/10.0
    z < 20.0 && return p10*exp(-lk_mf[]*(1.0 - exp(-α*(z - 10.0)))/α)
    p20 = p10*exp(-lk_mf[]*(1.0 - exp(-α*10.0))/α)
    return p20*exp(-lk_mf[]*(z - 20.0)/LK_TCOR)
end

@inline lk_density(z) = lk_mf[]*lk_pressure(z)/lk_temperature(z)

# Dipole field (B_x, B_z) at (x, z) for the signed strength B_d.
@inline function lk_dipole(x, z, Bd)
    s  = z - LK_ZD
    r2 = x*x + s*s
    r5 = r2*r2*sqrt(r2)
    return Bd*(3.0*s*s - r2)/r5, -3.0*Bd*s*x/r5
end

# Tube field (B_x, B_y, B_z) and its pressure deficit p₁ at (x, z).
@inline function lk_tube(x, z)
    ζ  = z - LK_ZT
    r2 = x*x + ζ*ζ
    q  = 1.0/LK_RT
    By = LK_BT*exp(-r2/LK_RT^2)
    p1 = 0.5*LK_BT^2*exp(-2.0*r2/LK_RT^2)*(q*q*(0.5*LK_RT^2 - r2) - 1.0)
    return q*ζ*By, By, -q*x*By, p1
end

# Φ_dip/Φ_tube of Paper II (Eqs. 9-10) on x = y = 0 for a reversed dipole of
# strength |B_d|: the total B_x changes sign at the separatrix z_sep above the axis.
# Φ_dip is the dipole flux above z_sep; Φ_tube the tube's from its axis to z_sep
# (def = "top", the printed Eq. 10) or the whole tube's azimuthal flux B_t R
# (def = "whole", the one Paper II's Fig. 1 follows). Returns (ratio, z_sep, Φ_tube, Φ_dip).
function lk_flux_ratio(Bdabs, ztop; n = 200001, def = "whole")
    z  = range(LK_ZT, ztop, length = n)
    Bx = [lk_tube(0.0, zz)[1] + lk_dipole(0.0, zz, -Bdabs)[1] for zz in z]
    i0 = findfirst(>(0.0), Bx)
    is = findnext(<(0.0), Bx, i0)
    is === nothing && return 0.0, ztop, 0.0, 0.0
    dz = step(z)
    Φt = def == "top" ? dz*(sum(Bx[1:is-1]) - 0.5*(Bx[1] + Bx[is-1])) : LK_BT*LK_RT
    Φd = dz*(sum(Bx[is:end]) - 0.5*(Bx[is] + Bx[end]))
    return abs(Φd)/Φt, z[is], Φt, Φd
end

# |B_d| giving the requested ratio (bisection; the ratio grows with |B_d|).
function lk_solve_dipole(ratio, ztop; def = "whole")
    ratio <= 0.0 && return 0.0
    lo, hi = 0.0, 1.0e5
    for _ = 1:60
        mid = 0.5*(lo + hi)
        lk_flux_ratio(mid, ztop; n = 20001, def = def)[1] < ratio ? (lo = mid) : (hi = mid)
    end
    return 0.5*(lo + hi)
end

function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && println(" Initialize fields for 3D GLM-MHD flux emergence (Leake et al. 2013, 2014) ........... ")

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "β", "Bmag", "log10ρ"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    inputs[:backend] != CPU()         && error(" problems/MHD/fluxEmergenceLeake2014: only the CPU backend is supported.")
    inputs[:SOL_VARS_TYPE] != TOTAL() && error(" problems/MHD/fluxEmergenceLeake2014: only SOL_VARS_TYPE = TOTAL() is supported.")

    lk_cz_a[]      = Float64(get(inputs, :lk_cz_gradient, 1.0))
    lk_mf[]        = Float64(get(inputs, :lk_mf, 1.0))
    lk_glm_alpha[] = Float64(get(inputs, :lk_glm_alpha, 1.0))
    lk_damp_xd[]   = Float64(get(inputs, :lk_damp_xd, 96.0))
    lk_damp_zd[]   = Float64(get(inputs, :lk_damp_zd, 180.0))
    lk_xmax[]      = MPI.Allreduce(maximum(abs, @view(mesh.coords[1,1:mesh.npoin])), MPI.MAX, comm)
    lk_zmax[]      = MPI.Allreduce(maximum(@view(mesh.coords[3,1:mesh.npoin])), MPI.MAX, comm)
    ltube          = get(inputs, :lk_tube, true)

    # Dipole: Paper II's flux ratio (default) or an explicit signed B_d.
    dip = String(get(inputs, :lk_dipole, "MD"))
    if haskey(inputs, :lk_Bd)
        lk_Bd[] = Float64(inputs[:lk_Bd])
    else
        haskey(LK_RATIO, dip) || error(" problems/MHD/fluxEmergenceLeake2014: :lk_dipole = \"$(dip)\", expected SD, MD, WD or ND.")
        fdef = String(get(inputs, :lk_flux_def, "whole"))
        fdef in ("whole", "top") || error(" problems/MHD/fluxEmergenceLeake2014: :lk_flux_def = \"$(fdef)\", expected whole or top.")
        lk_Bd[] = -lk_solve_dipole(LK_RATIO[dip], lk_zmax[]; def = fdef)     # reversed: B_x < 0 above the tube
    end
    ratio, zsep, Φt, Φd = lk_Bd[] < 0.0 ? lk_flux_ratio(-lk_Bd[], lk_zmax[]) : (0.0, NaN, NaN, NaN)
    rtop = lk_Bd[] < 0.0 ? lk_flux_ratio(-lk_Bd[], lk_zmax[]; def = "top")[1] : 0.0

    γm1 = γ_mhd - 1.0
    Bd  = lk_Bd[]
    ch_local = 0.0
    lk_seed  = parse(Float64, get(ENV, "LK_SEED", "0"))
    for ip = 1:mesh.npoin
        x, y, z = mesh.coords[1,ip], mesh.coords[2,ip], mesh.coords[3,ip]
        p0 = lk_pressure(z)
        ρ0 = lk_density(z)
        Bxd, Bzd = lk_dipole(x, z, Bd)

        # reference state: atmosphere + dipole
        q.qe[ip,1] = ρ0
        q.qe[ip,2] = 0.0; q.qe[ip,3] = 0.0; q.qe[ip,5] = 0.0
        q.qe[ip,4] = p0/γm1 + 0.5*(Bxd*Bxd + Bzd*Bzd)
        q.qe[ip,6] = Bxd; q.qe[ip,7] = 0.0; q.qe[ip,8] = Bzd
        q.qe[ip,9] = 0.0
        q.qe[ip,end] = p0

        Bxt, Byt, Bzt, p1 = ltube ? lk_tube(x, z) : (0.0, 0.0, 0.0, 0.0)
        p  = p0 + p1
        ρ  = ρ0*(1.0 + (p1/p0)*exp(-(y/LK_LAM)^2))
        Bx = Bxd + Bxt; By = Byt; Bz = Bzd + Bzt
        q.qn[ip,1] = ρ
        q.qn[ip,2] = 0.0; q.qn[ip,3] = 0.0; q.qn[ip,5] = 0.0
        q.qn[ip,4] = p/γm1 + 0.5*(Bx*Bx + By*By + Bz*Bz)
        q.qn[ip,6] = Bx; q.qn[ip,7] = By; q.qn[ip,8] = Bz
        q.qn[ip,9] = 0.0
        q.qn[ip,end] = p
        # LK_SEED=a (diagnostics): a vertical velocity of amplitude a, pseudo-random from node to
        # node (a hash of the coordinates, so that the copies of a shared node agree), in z ∈ [0, 30]
        if lk_seed > 0.0 && 0.0 <= z <= 30.0
            h = sin(12.9898*x + 78.233*y + 37.719*z)*43758.5453
            q.qn[ip,5]  = lk_seed*ρ*(2.0*(h - floor(h)) - 1.0)
            q.qn[ip,4] += 0.5*q.qn[ip,5]^2/ρ
        end

        ch_local = max(ch_local, sqrt((γ_mhd*p + Bx*Bx + By*By + Bz*Bz)/ρ))   # fast-speed bound, at rest
    end
    # GLM cleaning speed: the fastest wave of the initial state (the Alfvén speed
    # just above the transition region), held constant.
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    # Smallest nodal spacing (Dedner's damping rate α_p c_h/Δh).
    ngl = mesh.ngl
    dh_local = Inf
    for ie = 1:mesh.nelem
        a = mesh.connijk[ie,1,1,1]
        for b in (mesh.connijk[ie,2,1,1], mesh.connijk[ie,1,2,1], mesh.connijk[ie,1,1,2])
            dh_local = min(dh_local, sqrt(sum(abs2, @view(mesh.coords[:,b]) .- @view(mesh.coords[:,a]))))
        end
    end
    lk_glm_dh[] = MPI.Allreduce(dh_local, MPI.MIN, comm)

    get(inputs, :lwell_balanced, false) || rank == 0 &&
        println(" WARNING problems/MHD/fluxEmergenceLeake2014: :lwell_balanced is off; the discrete",
                " imbalance of the atmosphere and the dipole will drive flows from t = 0.")

    if rank == 0
        println(" Dipole $(dip): B_d = $(Bd) B₀L₀² (reversed), Φ_dip/Φ_tube = $(ratio) with the whole tube's flux $(Φt) ($(rtop) with its top half); z_sep = $(zsep), Φ_dip = $(Φd)")
        println(" Atmosphere (m_f = $(lk_mf[])): p(-30) = $(lk_pressure(-30.0)), ρ(-30) = $(lk_density(-30.0)); ρ(20) = $(lk_density(20.0)), ρ(top) = $(lk_density(lk_zmax[])); convection-zone gradient a = $(lk_cz_a[])")
        println(" Tube: B_t = $(LK_BT), R = $(LK_RT), q = 1/R, z_t = $(LK_ZT), λ = $(LK_LAM)$(ltube ? "" : " (OFF: :lk_tube => false)")")
        println(" GLM c_h = $(c_h_mhd[]); smallest nodal spacing Δh = $(lk_glm_dh[]); ψ damping rate α_p c_h/Δh = $(lk_glm_alpha[]*c_h_mhd[]/lk_glm_dh[]) per t₀")
        println(" Velocity damping: |x|, |y| > $(lk_damp_xd[]) (boundary at $(lk_xmax[])) and z > $(lk_damp_zd[]) (top at $(lk_zmax[]))")
        flush(stdout)
    end

    # Restart from a snapshot of this case (README, "Restart"). After c_h and Δh, which keep
    # the values of the t = 0 state, as in the run that wrote the snapshot.
    rdir = String(get(inputs, :lk_restart_dir, ""))
    isempty(rdir) || lk_restart_from_vtk!(q, mesh, inputs, rdir)

    rank == 0 && (println(" Initialize fields for 3D GLM-MHD flux emergence ........... DONE"); flush(stdout))
    return q
end

# Restart from the VTK output of an earlier run of this case: the snapshot holds the nine
# primitives (ρ, u, v, w, p, Bx, By, Bz, ψ) in Float64 at every mesh node, so the conserved
# state is rebuilt exactly, E = p/(γ−1) + ½ρ|v|² + ½|B|² + ½ψ². qe (atmosphere + dipole)
# and c_h are those of the analytic initial state. Each rank reads its own piece
# iter_N/iter_N_<rank+1>.vtu, so the restart needs the rank count of the run that wrote it.
# :lk_restart_time picks the snapshot nearest to it in <dir>/simulation.pvd (default: the
# last one), and inputs[:tinit] is set to its time. The DynSGS residual history starts
# afresh (its first residual is zero), so the continuation is not bitwise.
function lk_restart_from_vtk!(q, mesh, inputs, rdir::String)
    comm   = MPI.COMM_WORLD
    rank   = MPI.Comm_rank(comm)
    nparts = MPI.Comm_size(comm)
    pvd    = joinpath(rdir, "simulation.pvd")
    buf    = [NaN, -1.0]
    if rank == 0
        isfile(pvd) || error(" problems/MHD/fluxEmergenceLeake2014: restart: $(pvd) not found (:lk_restart_dir is the run's output directory).")
        ent = Tuple{Float64,Int}[]
        for line in readlines(pvd)
            m = match(r"timestep=\"([^\"]+)\"[^>]*file=\"iter_(\d+)\.pvtu\"", line)
            m === nothing || push!(ent, (parse(Float64, m[1]), parse(Int, m[2])))
        end
        isempty(ent) && error(" problems/MHD/fluxEmergenceLeake2014: restart: no snapshot listed in $(pvd).")
        want = get(inputs, :lk_restart_time, nothing)
        k    = want === nothing ? length(ent) : argmin([abs(e[1] - Float64(want)) for e in ent])
        buf .= (ent[k][1], Float64(ent[k][2]))
    end
    MPI.Bcast!(buf, 0, comm)
    tr, iout = buf[1], Int(buf[2])

    fname = joinpath(rdir, "iter_$(iout)", "iter_$(iout)_$(lpad(rank + 1, ndigits(nparts), '0')).vtu")
    isfile(fname) || error(" problems/MHD/fluxEmergenceLeake2014: restart: $(fname) not found; the restart needs the same number of MPI ranks as the run that wrote $(rdir).")
    names = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ"]
    vars, np = read_vtu_point_data(fname, vcat(names, ["__coords__"]))
    np == mesh.npoin || error(" problems/MHD/fluxEmergenceLeake2014: restart: $(fname) has $(np) points, this rank's mesh $(mesh.npoin): same mesh and rank count needed.")
    for nm in names
        haskey(vars, nm) || error(" problems/MHD/fluxEmergenceLeake2014: restart: field $(nm) missing in $(fname).")
    end
    ipmap = vtk_to_mesh_ipmap(vars["__coords__"], @view(mesh.coords[1,1:mesh.npoin]),
                              @view(mesh.coords[2,1:mesh.npoin]), @view(mesh.coords[3,1:mesh.npoin]))
    ρa = vars["ρ"]; ua = vars["u"]; va = vars["v"]; wa = vars["w"]; pa = vars["p"]
    Bxa = vars["Bx"]; Bya = vars["By"]; Bza = vars["Bz"]; ψa = vars["ψ"]
    γm1 = γ_mhd - 1.0
    nbad = 0
    for ip = 1:mesh.npoin
        j = ipmap[ip]
        ρ = ρa[j]; u = ua[j]; v = va[j]; w = wa[j]; p = pa[j]
        Bx = Bxa[j]; By = Bya[j]; Bz = Bza[j]; ψ = ψa[j]
        all(isfinite, (ρ, u, v, w, p, Bx, By, Bz, ψ)) || (nbad += 1)
        q.qn[ip,1] = ρ
        q.qn[ip,2] = ρ*u; q.qn[ip,3] = ρ*v; q.qn[ip,5] = ρ*w
        q.qn[ip,4] = p/γm1 + 0.5*ρ*(u*u + v*v + w*w) + 0.5*(Bx*Bx + By*By + Bz*Bz) + 0.5*ψ*ψ
        q.qn[ip,6] = Bx; q.qn[ip,7] = By; q.qn[ip,8] = Bz; q.qn[ip,9] = ψ
        q.qn[ip,end] = p
    end
    nbad = MPI.Allreduce(nbad, +, comm)
    nbad == 0 || error(" problems/MHD/fluxEmergenceLeake2014: restart: $(nbad) nodes of iter_$(iout) (t = $(tr)) are not finite.")
    inputs[:tinit] = tr
    if rank == 0
        println(" Restart: t = $(tr) from $(joinpath(rdir, "iter_$(iout)")) ($(nparts) piece$(nparts > 1 ? "s" : ""))")
        flush(stdout)
    end
    return nothing
end
