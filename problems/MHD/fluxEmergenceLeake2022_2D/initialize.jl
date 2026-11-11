#---------------------------------------------------------------------------------
# Initial state of the 2.5D analogue of the flux-emergence study of
#
#   J. E. Leake, M. G. Linton & S. K. Antiochos, "The role of reconnection in the
#   onset of solar eruptions", ApJ 934:10 (2022), Sec. 2.4
#
# in the plane perpendicular to the rope (README): the code's (x, y) are the paper's
# (x, z), y up; the rope axis, the paper's y, is the invariant (out-of-plane, "z")
# direction. Units of user_flux.jl (L₀, B₀, g₀; p = ρT/m_f with m_f = 1, g = 1,
# magnetic pressure ½|B|²).
#
# 1. Atmosphere (the paper's Eqs. 23-25, as problems/MHD/fluxEmergenceLeake2022),
#    hydrostatic, ρ(0) = 1, with h the height (the code's y):
#      T = 1 − (γ−1)/γ · h      h ≤ 0     (convection zone, polytropic)
#      T = 1                    0 < h < 10 (photosphere/chromosphere)
#      T = 150^((h−10)/10)      10 ≤ h < 20 (transition region)
#      T = 150                  h ≥ 20    (corona)
#    dp/dh = −ρ = −m_f p/T in closed form per layer (lk2d_pressure).
#
# 2. Coronal field: a LINE dipole at (0, h_d) below the domain, invariant along the
#    rope, the 2D potential field
#      A = −D s/r²,  B = (∂A/∂y, −∂A/∂x) = D ((s² − x²)/r⁴, −2 s x/r⁴, 0),
#      s = h − h_d, r² = x² + s²,
#    pointing along +x above the dipole, as the paper's point dipole (moment −x̂)
#    does above the rope. h_d = −100 and D = 150 give |B| = D/h_d² = 0.015 B₀ =
#    19.5 G at (0, 0), the 3D dipole's value there (the 3D field falls off as 1/s³,
#    this one as 1/s²: README). The paper's point dipole (3D, falling as 1/r³ in
#    every direction) has no invariant counterpart, and A ∝ s/r³ (the 2014 cases'
#    formula) is not potential in 2D (∇²A = 3s/r⁵ ≠ 0, a current in the corona).
#
# 3. Flux rope (the paper's Eqs. 28-34) along the out-of-plane direction at
#    (0, h_rope = −12), width w = 2.5, strength B_t = 5, right-handed twist:
#      B_axial = σ B_t e^(−r²/w²)   (the code's Bz),
#      B_twist = (r/w) B_axial      (in the plane): Bx = σ ζ B_a/w, By = −σ x B_a/w,
#    r² = x² + ζ², ζ = h − h_rope, B_a = B_t e^(−r²/w²). σ = +1 is the paper's θ = 0
#    (the twist above the axis parallel to the dipole field: stable), σ = −1 its
#    θ = π, simulation 4 (antiparallel: erupts in 3D); the paper's other angles
#    have no invariant counterpart. Pressure deficit (Lorentz-force balance)
#      p₁ = ½ B_t² e^(−2r²/w²) [(w²/2 − r²)/w² − 1],
#    p = p₀(h) + p₁ and ρ = ρ₀(h) [1 + b p₁/p₀(h)]: b = 1 (:lk2d_buoyancy) is the
#    middle of the 3D rope (y_r = 0, isothermal with its surroundings), buoyant
#    along its whole length here; b = 0 the neutrally buoyant ends.
#
# The reference state qe (well balance, DynSGS departure, positivity floors) is the
# atmosphere plus the dipole, without the rope.
#---------------------------------------------------------------------------------

if !@isdefined(LK2D_TCOR)
    const LK2D_TCOR = 150.0    # T_cor/T_ph
    const LK2D_HR   = -12.0    # rope axis height h_rope
    const LK2D_W    = 2.5      # rope width w
    const LK2D_BT   = 5.0      # rope axial field B_t
end
if !@isdefined(lk2d_D)
    const lk2d_D         = Ref{Float64}(150.0)    # line-dipole strength (B_x = +D/s² above it)
    const lk2d_hd        = Ref{Float64}(-100.0)   # line-dipole depth h_d
    const lk2d_σ         = Ref{Float64}(-1.0)     # rope orientation: +1 parallel (θ = 0), −1 antiparallel (θ = π)
    const lk2d_b         = Ref{Float64}(1.0)      # buoyancy factor b
    const lk2d_cz_a      = Ref{Float64}(1.0)      # convection-zone gradient factor (1: Eq. 23 as printed)
    const lk2d_mf        = Ref{Float64}(1.0)      # m_f of p = ρT/m_f
    const lk2d_glm_alpha = Ref{Float64}(1.0)      # Dedner damping α_p (rate α_p c_h/Δh)
    const lk2d_glm_dh    = Ref{Float64}(1.0)
    const lk2d_damp_xd   = Ref{Float64}(207.9)    # side damping from |x| = 3/8 of the width
    const lk2d_damp_yd   = Ref{Float64}(464.2)    # top damping from y = 7/8 of the top height
    const lk2d_xmax      = Ref{Float64}(277.2)
    const lk2d_ymax      = Ref{Float64}(530.5)
end

# Temperature and pressure of the background atmosphere at height h (item 1).
@inline function lk2d_temperature(h)
    a = lk2d_cz_a[]
    h <= 0.0  && return 1.0 - a*(γ_mhd - 1.0)/γ_mhd*h
    h < 10.0  && return 1.0
    h < 20.0  && return LK2D_TCOR^((h - 10.0)/10.0)
    return LK2D_TCOR
end

@inline function lk2d_pressure(h)
    a  = lk2d_cz_a[]
    p0 = 1.0/lk2d_mf[]                                # ρ(0) = 1, T(0) = 1
    if h <= 0.0
        return p0*lk2d_temperature(h)^(lk2d_mf[]*γ_mhd/(a*(γ_mhd - 1.0)))
    end
    p10 = p0*exp(-lk2d_mf[]*10.0)
    h < 10.0 && return p0*exp(-lk2d_mf[]*h)
    α = log(LK2D_TCOR)/10.0
    h < 20.0 && return p10*exp(-lk2d_mf[]*(1.0 - exp(-α*(h - 10.0)))/α)
    p20 = p10*exp(-lk2d_mf[]*(1.0 - exp(-α*10.0))/α)
    return p20*exp(-lk2d_mf[]*(h - 20.0)/LK2D_TCOR)
end

@inline lk2d_density(h) = lk2d_mf[]*lk2d_pressure(h)/lk2d_temperature(h)

# Field of the line dipole (item 2): (Bx, By) at (x, h).
@inline function lk2d_dipole(x, h)
    D  = lk2d_D[]
    s  = h - lk2d_hd[]
    r2 = x*x + s*s
    r4 = r2*r2
    return D*(s*s - x*x)/r4, -2.0*D*s*x/r4
end

# Rope field (Bx, By, Bz) and pressure deficit p₁ at (x, h) (item 3).
@inline function lk2d_rope(x, h)
    σ  = lk2d_σ[]
    ζ  = h - LK2D_HR
    r2 = x*x + ζ*ζ
    w2 = LK2D_W*LK2D_W
    Ba = LK2D_BT*exp(-r2/w2)
    p1 = 0.5*LK2D_BT^2*exp(-2.0*r2/w2)*((0.5*w2 - r2)/w2 - 1.0)
    return σ*ζ*Ba/LK2D_W, -σ*x*Ba/LK2D_W, σ*Ba, p1
end

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && println(" Initialize fields for 2.5D GLM-MHD flux emergence (analogue of Leake et al. 2022) ........... ")

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "β", "Bmag"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    inputs[:backend] != CPU()         && error(" problems/MHD/fluxEmergenceLeake2022_2D: only the CPU backend is supported.")
    inputs[:SOL_VARS_TYPE] != TOTAL() && error(" problems/MHD/fluxEmergenceLeake2022_2D: only SOL_VARS_TYPE = TOTAL() is supported.")

    orient = String(get(inputs, :lk2d_orientation, "antiparallel"))
    orient in ("antiparallel", "parallel") || error(" problems/MHD/fluxEmergenceLeake2022_2D: :lk2d_orientation = \"$(orient)\", expected antiparallel or parallel.")
    lk2d_σ[]         = orient == "parallel" ? 1.0 : -1.0
    ldip             = get(inputs, :lk2d_dipole, true)
    lk2d_D[]         = ldip ? Float64(get(inputs, :lk2d_D, 150.0)) : 0.0
    lk2d_hd[]        = Float64(get(inputs, :lk2d_hd, -100.0))
    lk2d_b[]         = Float64(get(inputs, :lk2d_buoyancy, 1.0))
    lk2d_cz_a[]      = Float64(get(inputs, :lk2d_cz_gradient, 1.0))
    lk2d_mf[]        = Float64(get(inputs, :lk2d_mf, 1.0))
    lk2d_glm_alpha[] = Float64(get(inputs, :lk2d_glm_alpha, 1.0))
    lk2d_xmax[]      = MPI.Allreduce(maximum(abs, @view(mesh.coords[1,1:mesh.npoin]); init = 0.0), MPI.MAX, comm)
    ymin             = MPI.Allreduce(minimum(@view(mesh.coords[2,1:mesh.npoin]); init = Inf), MPI.MIN, comm)
    lk2d_ymax[]      = MPI.Allreduce(maximum(@view(mesh.coords[2,1:mesh.npoin]); init = -Inf), MPI.MAX, comm)
    # damping regions (the paper's Eqs. 21-22): |x| > 3L_x/8 (L_x the full width, i.e. the
    # outer eighth of each side) and y > 7 max(y)/8
    lk2d_damp_xd[]   = Float64(get(inputs, :lk2d_damp_x_frac, 0.75))*lk2d_xmax[]
    lk2d_damp_yd[]   = Float64(get(inputs, :lk2d_damp_z_frac, 7.0/8.0))*lk2d_ymax[]
    lk2d_hd[] >= ymin && error(" problems/MHD/fluxEmergenceLeake2022_2D: the dipole (h_d = $(lk2d_hd[])) must lie below the domain (y ≥ $(ymin)).")
    lrope = get(inputs, :lk2d_rope, true)

    γm1 = γ_mhd - 1.0
    ch_local = 0.0
    for ip = 1:mesh.npoin
        x, h = mesh.coords[1,ip], mesh.coords[2,ip]
        p0 = lk2d_pressure(h)
        ρ0 = lk2d_density(h)
        Bxd, Byd = lk2d_dipole(x, h)

        # reference state: atmosphere + dipole
        q.qe[ip,1] = ρ0
        q.qe[ip,2] = 0.0; q.qe[ip,3] = 0.0; q.qe[ip,5] = 0.0
        q.qe[ip,4] = p0/γm1 + 0.5*(Bxd*Bxd + Byd*Byd)
        q.qe[ip,6] = Bxd; q.qe[ip,7] = Byd; q.qe[ip,8] = 0.0
        q.qe[ip,9] = 0.0
        q.qe[ip,end] = p0

        Bxr, Byr, Bzr, p1 = lrope ? lk2d_rope(x, h) : (0.0, 0.0, 0.0, 0.0)
        p  = p0 + p1
        ρ  = ρ0*(1.0 + lk2d_b[]*p1/p0)
        Bx = Bxd + Bxr; By = Byd + Byr; Bz = Bzr
        q.qn[ip,1] = ρ
        q.qn[ip,2] = 0.0; q.qn[ip,3] = 0.0; q.qn[ip,5] = 0.0
        q.qn[ip,4] = p/γm1 + 0.5*(Bx*Bx + By*By + Bz*Bz)
        q.qn[ip,6] = Bx; q.qn[ip,7] = By; q.qn[ip,8] = Bz
        q.qn[ip,9] = 0.0
        q.qn[ip,end] = p

        ch_local = max(ch_local, sqrt((γ_mhd*p + Bx*Bx + By*By + Bz*Bz)/ρ))   # fast-speed bound, at rest
    end
    # GLM cleaning speed: the fastest wave of the initial state (the Alfvén speed just
    # above the transition region), held constant.
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    # Smallest nodal spacing (Dedner's damping rate α_p c_h/Δh).
    dh_local = Inf
    for ie = 1:mesh.nelem
        a = mesh.connijk[ie,1,1]
        for b in (mesh.connijk[ie,2,1], mesh.connijk[ie,1,2])
            dh_local = min(dh_local, sqrt((mesh.coords[1,b] - mesh.coords[1,a])^2 + (mesh.coords[2,b] - mesh.coords[2,a])^2))
        end
    end
    lk2d_glm_dh[] = MPI.Allreduce(dh_local, MPI.MIN, comm)

    get(inputs, :lwell_balanced, false) || rank == 0 &&
        println(" WARNING problems/MHD/fluxEmergenceLeake2022_2D: :lwell_balanced is off; the discrete",
                " imbalance of the atmosphere and the dipole will drive flows from t = 0.")

    if rank == 0
        println(" Rope: $(orient) (σ = $(lk2d_σ[]); the paper's θ = $(orient == "parallel" ? "0" : "π")), B_t = $(LK2D_BT), w = $(LK2D_W), h_rope = $(LK2D_HR), buoyancy b = $(lk2d_b[])$(lrope ? "" : " (OFF: :lk2d_rope => false)")")
        println(" Line dipole: $(ldip ? "D = $(lk2d_D[]) at h_d = $(lk2d_hd[]); B_x at (0, 0) = $(lk2d_D[]/lk2d_hd[]^2) B₀ = $(round(1300*lk2d_D[]/lk2d_hd[]^2, digits=2)) G" : "OFF (:lk2d_dipole => false)")")
        println(" Atmosphere (m_f = $(lk2d_mf[])): p(-30) = $(lk2d_pressure(-30.0)), ρ(-30) = $(lk2d_density(-30.0)); ρ(20) = $(lk2d_density(20.0)), ρ(top) = $(lk2d_density(lk2d_ymax[]))")
        println(" GLM c_h = $(c_h_mhd[]); smallest nodal spacing Δh = $(lk2d_glm_dh[]); ψ damping rate α_p c_h/Δh = $(lk2d_glm_alpha[]*c_h_mhd[]/lk2d_glm_dh[]) per t₀")
        println(" Velocity damping: |x| > $(lk2d_damp_xd[]) (boundary at $(lk2d_xmax[])) and y > $(lk2d_damp_yd[]) (top at $(lk2d_ymax[]))")
        flush(stdout)
    end

    # Restart from a snapshot of this case (README, "Restart"). After c_h and Δh, which keep
    # the values of the t = 0 state, as in the run that wrote the snapshot.
    rdir = String(get(inputs, :lk2d_restart_dir, ""))
    isempty(rdir) || lk2d_restart_from_vtk!(q, mesh, inputs, rdir)

    rank == 0 && (println(" Initialize fields for 2.5D GLM-MHD flux emergence ........... DONE"); flush(stdout))
    return q
end

# Restart from the VTK output of an earlier run of this case (as the 3D Leake cases): the
# snapshot holds the nine primitives (ρ, u, v, w, p, Bx, By, Bz, ψ) in Float64 at every mesh
# node, so the conserved state is rebuilt exactly, E = p/(γ−1) + ½ρ|v|² + ½|B|² + ½ψ². qe and
# c_h are those of the analytic initial state. Each rank reads its own piece
# iter_N/iter_N_<rank+1>.vtu, so the restart needs the rank count of the run that wrote it.
# :lk2d_restart_time picks the snapshot nearest to it in <dir>/simulation.pvd (default: the
# last one), and inputs[:tinit] is set to its time. The DynSGS residual history starts
# afresh (its first residual is zero), so the continuation is not bitwise.
function lk2d_restart_from_vtk!(q, mesh, inputs, rdir::String)
    comm   = MPI.COMM_WORLD
    rank   = MPI.Comm_rank(comm)
    nparts = MPI.Comm_size(comm)
    pvd    = joinpath(rdir, "simulation.pvd")
    buf    = [NaN, -1.0]
    if rank == 0
        isfile(pvd) || error(" problems/MHD/fluxEmergenceLeake2022_2D: restart: $(pvd) not found (:lk2d_restart_dir is the run's output directory).")
        ent = Tuple{Float64,Int}[]
        for line in readlines(pvd)
            m = match(r"timestep=\"([^\"]+)\"[^>]*file=\"iter_(\d+)\.pvtu\"", line)
            m === nothing || push!(ent, (parse(Float64, m[1]), parse(Int, m[2])))
        end
        isempty(ent) && error(" problems/MHD/fluxEmergenceLeake2022_2D: restart: no snapshot listed in $(pvd).")
        want = get(inputs, :lk2d_restart_time, nothing)
        k    = want === nothing ? length(ent) : argmin([abs(e[1] - Float64(want)) for e in ent])
        buf .= (ent[k][1], Float64(ent[k][2]))
    end
    MPI.Bcast!(buf, 0, comm)
    tr, iout = buf[1], Int(buf[2])

    fname = joinpath(rdir, "iter_$(iout)", "iter_$(iout)_$(lpad(rank + 1, ndigits(nparts), '0')).vtu")
    isfile(fname) || error(" problems/MHD/fluxEmergenceLeake2022_2D: restart: $(fname) not found; the restart needs the same number of MPI ranks as the run that wrote $(rdir).")
    names = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ"]
    vars, np = read_vtu_point_data(fname, vcat(names, ["__coords__"]))
    np == mesh.npoin || error(" problems/MHD/fluxEmergenceLeake2022_2D: restart: $(fname) has $(np) points, this rank's mesh $(mesh.npoin): same mesh and rank count needed.")
    for nm in names
        haskey(vars, nm) || error(" problems/MHD/fluxEmergenceLeake2022_2D: restart: field $(nm) missing in $(fname).")
    end
    # the 2D output writes its points with z = 0
    ipmap = vtk_to_mesh_ipmap(vars["__coords__"], @view(mesh.coords[1,1:mesh.npoin]),
                              @view(mesh.coords[2,1:mesh.npoin]), zeros(mesh.npoin))
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
    nbad == 0 || error(" problems/MHD/fluxEmergenceLeake2022_2D: restart: $(nbad) nodes of iter_$(iout) (t = $(tr)) are not finite.")
    inputs[:tinit] = tr
    if rank == 0
        println(" Restart: t = $(tr) from $(joinpath(rdir, "iter_$(iout)")) ($(nparts) piece$(nparts > 1 ? "s" : ""))")
        flush(stdout)
    end
    return nothing
end
