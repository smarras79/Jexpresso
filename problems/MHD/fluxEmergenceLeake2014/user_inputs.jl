#---------------------------------------------------------------------------------
# 3D emergence of a twisted convection-zone flux tube into a coronal dipole:
# Leake, Linton & Antiochos, ApJ 787:46 (2014), simulations SD, MD, WD, with the
# setup of Leake, Linton & Török, ApJ 778:99 (2013). See README.md.
#
#   julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2014
#   mpiexec -n 128 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2014
#
# LK_REF picks the mesh (LK.geo) and the time step:
#   1: 16³ elements  (smoke test; 19 of the paper's cells per element)
#   2: 38³ elements  (half the paper's resolution at :nop => 4)
#   3: 76³ elements  (the paper's 304³ cells at :nop => 4; generate it first:
#                     gmsh -3 LK.geo -setnumber ref 3 -o LK_76x76x76.msh)
#---------------------------------------------------------------------------------

# Positivity floors as the Runge-Kutta stage limiter (as problems/MHD/fluxEmergenceSon2025):
# ρ ≥ f ρ_e(z) and p ≥ f p_e(z) of the reference atmosphere, p raised through E at fixed
# velocity and field. Not conservative where it fires; the hits are counted over all ranks and
# reported by rank 0 each time the total passes a power of ten.
if !@isdefined(lk_floor_frac)
    const lk_floor_frac = Ref{Float64}(0.2)
    const lk_floor_hits = Ref{Int}(0)
    const lk_floor_next = Ref{Int}(1)
end
# CFL monitor (the time step is fixed, the wave speeds are not): every 100 stages the largest
# (|v| + c_f) Δt/Δh over the nodes, with c_f = √((γp + |B|²)/ρ) and Δh the smallest LGL gap of
# the elements around the node, is reported by rank 0 when it has grown by 20 %, as a
# WARNING above 0.5 (the step is set for ≈ 0.25 initially: README).
if !@isdefined(lk_cfl_dh)
    const lk_cfl_dh    = Ref{Vector{Float64}}(Float64[])
    const lk_cfl_count = Ref{Int}(0)
    const lk_cfl_shown = Ref{Float64}(0.0)
end

function lk_node_spacing(mesh)
    ngl = Int(mesh.ngl)
    h   = fill(Inf, Int(mesh.npoin))
    for ie = 1:Int(mesh.nelem)
        he = Inf
        for m = 1:ngl-1, (a, b) in ((mesh.connijk[ie,m,1,1], mesh.connijk[ie,m+1,1,1]),
                                    (mesh.connijk[ie,1,m,1], mesh.connijk[ie,1,m+1,1]),
                                    (mesh.connijk[ie,1,1,m], mesh.connijk[ie,1,1,m+1]))
            he = min(he, sqrt((mesh.coords[1,b] - mesh.coords[1,a])^2 + (mesh.coords[2,b] - mesh.coords[2,a])^2 +
                              (mesh.coords[3,b] - mesh.coords[3,a])^2))
        end
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[ie,i,j,k]
            h[ip] = min(h[ip], he)
        end
    end
    return h
end

function lk_positivity_limiter!(u, integrator, p, t)
    npoin = p.mesh.npoin
    qe    = p.qp.qe
    frac  = lk_floor_frac[]
    γm1   = γ_mhd - 1.0
    nhit  = 0
    lk_cfl_count[] += 1
    docfl = lk_cfl_count[] % 100 == 1
    if docfl && length(lk_cfl_dh[]) != npoin
        lk_cfl_dh[] = lk_node_spacing(p.mesh)
    end
    dh    = lk_cfl_dh[]
    Δt    = Float64(p.Δt)
    cfl = 0.0; vmx = 0.0; cmx = 0.0
    @inbounds for ip = 1:npoin
        iρ = ip; iρu = npoin + ip; iρv = 2npoin + ip; iE = 3npoin + ip; iρw = 4npoin + ip
        iBx = 5npoin + ip; iBy = 6npoin + ip; iBz = 7npoin + ip; iψ = 8npoin + ip
        ρmin = frac*qe[ip,1]
        pmin = frac*qe[ip,end]
        ρ = u[iρ]
        if ρ < ρmin
            u[iρ] = ρmin; ρ = ρmin; nhit += 1
        end
        ke = 0.5*(u[iρu]^2 + u[iρv]^2 + u[iρw]^2)/ρ
        me = 0.5*(u[iBx]^2 + u[iBy]^2 + u[iBz]^2) + 0.5*u[iψ]^2
        if γm1*(u[iE] - ke - me) < pmin
            u[iE] = pmin/γm1 + ke + me; nhit += 1
        end
        if docfl
            B2 = u[iBx]^2 + u[iBy]^2 + u[iBz]^2
            v  = sqrt(u[iρu]^2 + u[iρv]^2 + u[iρw]^2)/ρ
            cf = sqrt(max(γ_mhd*γm1*(u[iE] - ke - me) + B2, 0.0)/ρ)
            c  = (v + cf)*Δt/dh[ip]
            if c > cfl
                cfl = c; vmx = v; cmx = cf
            end
        end
    end
    if docfl
        comm = MPI.COMM_WORLD
        g = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(cfl, MPI.MAX, comm) : cfl
        if g > 1.2*lk_cfl_shown[]
            if MPI.Comm_rank(comm) == 0       # the |v|, c_f printed are those of rank 0's largest
                print(string(g > 0.5 ? " # WARNING CFL: " : " # CFL: ", "max (|v| + c_f) Δt/Δh = ", round(g, sigdigits = 3),
                             " at t = ", round(t, digits = 3), " (rank 0: |v| = ", round(vmx, sigdigits = 3), ", c_f = ",
                             round(cmx, sigdigits = 3), ")", g > 0.5 ? "; reduce :Δt" : "", "\n"))
                flush(stdout)
            end
            lk_cfl_shown[] = g
        end
    end
    # One report for all ranks, from rank 0, when the total passes 1, 10, 100, ...: every rank
    # calls the stage limiter at every Runge-Kutta stage, so the sum is a safe collective. t is
    # the stage time t_n + c_i Δt (in t₀), not an output time. Nodes on a partition interface
    # are counted once per rank that holds them.
    comm = MPI.COMM_WORLD
    nall = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(nhit, +, comm) : nhit
    lk_floor_hits[] += nall
    if nall > 0 && lk_floor_hits[] >= lk_floor_next[]
        if MPI.Comm_rank(comm) == 0
            print(string(" # positivity floors: ", lk_floor_hits[], " node hits so far on all ranks (",
                         nall, " in this stage) at t = ", round(t, digits = 3), "\n"))
            flush(stdout)
        end
        while lk_floor_next[] <= lk_floor_hits[]
            lk_floor_next[] *= 10
        end
    end
    lk_diag!(u, p, t)
    return nothing
end

# TEMP DIAGNOSTIC (LK_DIAG=n env var: every n-th stage): global maxima of the state.
if !@isdefined(lk_diag_calls)
    const lk_diag_calls = Ref{Int}(0)
end
function lk_diag!(u, p, t)
    n = parse(Int, get(ENV, "LK_DIAG", "0"))
    n <= 0 && return nothing
    lk_diag_calls[] += 1
    lk_diag_calls[] % n == 0 || return nothing
    npoin = p.mesh.npoin; qe = p.qp.qe; γm1 = γ_mhd - 1.0; X = p.mesh.coords
    vmax = 0.0; wmax = 0.0; ρrmin = Inf; prmax = 0.0; prmin = Inf; ψmax = 0.0; Tmax = 0.0; dBmax = 0.0
    nρfl = 0; npfl = 0; zfl_lo = Inf; zfl_hi = -Inf
    xw = (0.0, 0.0, 0.0); xv = xw; xρ = xw; xp = xw; xpm = xw
    @inbounds for ip = 1:npoin
        ρ = u[ip]; ru = u[npoin+ip]; rv = u[2npoin+ip]; E = u[3npoin+ip]; rw = u[4npoin+ip]
        Bx = u[5npoin+ip]; By = u[6npoin+ip]; Bz = u[7npoin+ip]; ψ = u[8npoin+ip]
        B2 = Bx*Bx + By*By + Bz*Bz
        v = sqrt(ru*ru + rv*rv + rw*rw)/ρ
        pr = γm1*(E - 0.5*(ru*ru + rv*rv + rw*rw)/ρ - 0.5*B2 - 0.5*ψ*ψ)
        cf = sqrt(max((γ_mhd*pr + B2)/ρ, 0.0))
        v > vmax && (vmax = v; xv = (X[1,ip], X[2,ip], X[3,ip]))
        v + cf > wmax && (wmax = v + cf; xw = (X[1,ip], X[2,ip], X[3,ip]))
        rr = ρ/qe[ip,1]; pp = pr/qe[ip,end]
        rr < ρrmin && (ρrmin = rr; xρ = (X[1,ip], X[2,ip], X[3,ip]))
        pp > prmax && (prmax = pp; xp = (X[1,ip], X[2,ip], X[3,ip]))
        pp < prmin && (prmin = pp; xpm = (X[1,ip], X[2,ip], X[3,ip]))
        rr <= 1.000001*lk_floor_frac[] && (nρfl += 1)
        if rr <= 1.000001*lk_floor_frac[] || pp <= 1.000001*lk_floor_frac[]
            pp <= 1.000001*lk_floor_frac[] && (npfl += 1)
            zfl_lo = min(zfl_lo, X[3,ip]); zfl_hi = max(zfl_hi, X[3,ip])
        end
        ψmax = max(ψmax, abs(ψ)); Tmax = max(Tmax, lk_mf[]*pr/ρ)
        dBmax = max(dBmax, sqrt((Bx-qe[ip,6])^2 + By^2 + (Bz-qe[ip,8])^2))
    end
    # DynSGS + explicit coefficient (kinematic) of the momentum and field slots, and the
    # explicit-diffusion number μ Δt/Δh² with Δh the element's smallest LGL gap (the
    # Laplacian's largest eigenvalue is ≈ 5.5 μ/Δh² per direction at N = 4 and the
    # real-axis limit of CarpenterKennedy2N54 is 4.66: trouble from D ≈ 0.3 on).
    μd = p.μ_dsgs; Δt = Float64(p.Δt); dh = lk_cfl_dh[]
    μρ = 0.0; μm = 0.0; μB = 0.0; Dmax = 0.0; xD = (0.0, 0.0, 0.0)
    if size(μd, 1) == p.mesh.nelem && length(dh) == npoin
        @inbounds for ie = 1:p.mesh.nelem
            μρ = max(μρ, μd[ie,1]); μm = max(μm, μd[ie,2]); μB = max(μB, μd[ie,6])
            ip = p.mesh.connijk[ie,2,2,2]
            D = maximum(@view μd[ie,1:8])*Δt/dh[ip]^2
            D > Dmax && (Dmax = D; xD = (X[1,ip], X[2,ip], X[3,ip]))
        end
    end
    comm = MPI.COMM_WORLD
    loc = [vmax, wmax, -ρrmin, prmax, -prmin, ψmax, Tmax, dBmax, μρ, μm, μB, Dmax, -zfl_lo, zfl_hi]
    glo = MPI.Allreduce(loc, MPI.MAX, comm)
    cnt = MPI.Allreduce([nρfl, npfl], +, comm)
    if MPI.Comm_rank(comm) == 0
        println(" #DIAG t=", round(t, digits=4), " |v|max=", round(glo[1], sigdigits=4),
                " (|v|+cf)max=", round(glo[2], sigdigits=4), " ρ/ρe min=", round(-glo[3], sigdigits=3),
                " p/pe max=", round(glo[4], sigdigits=4), " p/pe min=", round(-glo[5], sigdigits=3),
                " |ψ|max=", round(glo[6], sigdigits=3), " Tmax=", round(glo[7], sigdigits=4),
                " |δB|max=", round(glo[8], sigdigits=3), " nodes at ρ/p floor=", cnt[1], "/", cnt[2],
                " μmax(ρ,m,B)=", round(glo[9], sigdigits=3), ",", round(glo[10], sigdigits=3), ",",
                round(glo[11], sigdigits=3), " μΔt/Δh²max=", round(glo[12], sigdigits=3),
                cnt[1] + cnt[2] > 0 ? string(" floors at z=", round(-glo[13], digits=1), "..", round(glo[14], digits=1)) : "")
    end
    # where: each extreme is printed by the rank that holds it
    for (val, g, x, name) in ((wmax, glo[2], xw, "(|v|+cf)max"), (vmax, glo[1], xv, "|v|max"),
                              (-ρrmin, glo[3], xρ, "ρ/ρe min"), (prmax, glo[4], xp, "p/pe max"),
                              (-prmin, glo[5], xpm, "p/pe min"),
                              (Dmax, glo[12], xD, "μΔt/Δh² max"))
        val == g && println(" #DIAG   ", name, " at ", round.(x, digits=2), " (rank ", MPI.Comm_rank(comm), ")")
    end
    flush(stdout)
    return nothing
end

# Explicit dissipation of the papers (Re = Rm = 100, i.e. 0.01 in units of L₀v₀),
# per element centroid, added to the DynSGS coefficient of the momentum (ν) and
# field (η) slots by compute_dsgs_viscosity!(::DSGS_MHD, ::NSD_3D):
#   η = 0.01 [tanh(−(r_η − 100)/5) + 1]/2,  r_η = √(x² + y²)   (Paper I Eq. 9)
#   ν = 0.01, kinematic (README: the papers' dynamic ν would give Re ~ 10⁻⁶ in the corona)
lk_dissipation(x, y, z) = (0.01, 0.01*0.5*(tanh(-(sqrt(x*x + y*y) - 100.0)/5.0) + 1.0))

function user_inputs()

    LK_REF = parse(Int, get(ENV, "LK_REF", "2"))
    mesh   = ("LK_16x16x16.msh", "LK_38x38x38.msh", "LK_76x76x76.msh")[LK_REF]
    # Δt for (|v| + c_f) Δt/Δh ≈ 0.25 at |v| + c_f ≈ 120 on the smallest LGL gap Δh of each mesh:
    # the fast speed just above the transition region is 92 (SD), 72 (MD), 50 (WD) at t = 0, and
    # Paper II's reconnection outflows reach 30 (README). The CFL monitor reports the actual value.
    Δt     = (3.0e-3, 1.2e-3, 6.0e-4)[LK_REF]
    # Environment overrides (README, "Restart"): LK_TEND sets :tend, LK_OUT_DT a snapshot every
    # LK_OUT_DT from 0 to :tend (only those after the start are written), LK_RESTART_DIR
    # (+ LK_RESTART_TIME) restarts from a snapshot of an earlier run.
    tend   = parse(Float64, get(ENV, "LK_TEND", "160.0"))
    dtout  = parse(Float64, get(ENV, "LK_OUT_DT", "5.0"))

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(stage_limiter! = lk_positivity_limiter!,
                                                      williamson_condition = false),
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => tend,       # 160: Paper II Figs. 10-11, MD erupts at 120-160 t₀
        :diagnostics_at_times => (0.0:dtout:tend),
        :conservation_every   => 10,
        :conservation_slots   => [1, 4],
        :lstep_heartbeat      => true,       # "#   step N   t = …" every 100 steps: snapshots are 1667 steps apart at LK_REF = 1
        :restart_time         => 0.0,
        :lrestart             => false,
        :lk_restart_dir       => get(ENV, "LK_RESTART_DIR", ""),   # "" = start from the initial state
        :lsource              => true,       # gravity, ψ damping, velocity damping
        :lwell_balanced       => true,       # subtract the discrete RHS of qe (atmosphere + dipole)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # The papers' problem (README):
        #---------------------------------------------------------------------------
        :lk_dipole            => "MD",       # SD, MD, WD (Φ_dip/Φ_tube = 0.13, 0.1, 0.067) or ND
        :lk_flux_def          => "whole",    # Φ_tube of the whole tube (Paper II's Fig. 1; README) or "top" (its Eq. 10)
        :lk_tube              => true,
        :lk_cz_gradient       => 1.0,        # a of Paper I Eq. 11 (1 = exactly adiabatic with m_f = 1)
        :lk_mf                => 1.0,        # p = ρT/m_f: 1 reproduces the papers' Fig. 1 (README); 1.25 = literal text
        :lk_glm_alpha         => 1.0,        # ψ damping rate α_p c_h/Δh (README: GLM transient)
        :lk_damp_xd           => 96.0,
        :lk_damp_zd           => 180.0,      # Paper II (Paper I: 170)
        :mhd_dissipation      => lk_dissipation,
        #---------------------------------------------------------------------------
        # DynSGS in conserved form on the departure from the reference state, as the
        # 2D flux-emergence case: legacy sensor, element norms floored at the local
        # scales, fast-speed floors (β ~ 10⁻² in the corona), thermal-energy residual.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model           => DSGS_MHD(),
        :dsgs_sensor          => "legacy",
        :dsgs_CR              => 1.0,
        :dsgs_Cmax            => 0.5,
        :dsgs_Cmin            => 0.0,
        :dsgs_gamma           => 5.0/3.0,
        :dsgs_Prt             => 1.0,
        :dsgs_norms           => "element",
        :dsgs_local_rel       => 1.0,
        :dsgs_conserved       => true,
        # Reference-weighted operator (user_primitives.jl, README): rho, rho v, E diffused as
        # (q - q_e)/rho_e with mu rho_e, which keeps the transition region from being evacuated and
        # heated by the diffused departure of the dense side (LK_REF_WEIGHT=0: the absolute form).
        :dsgs_ref_weight      => get(ENV, "LK_REF_WEIGHT", "0") == "1",   # being tested in 3D: LK_REF_WEIGHT=1
        :dsgs_fast_floors     => true,
        :dsgs_thermal_residual => true,
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh: LK.geo (the papers' stretching, element boundaries every 4 cells at ref 3)
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        # Balanced 3D (METIS) partition. The default column partition bins the x-y bounding box
        # uniformly, which on this horizontally stretched mesh leaves ranks without elements from
        # a few hundred ranks on and overloads the centre at any count (README, Notes).
        :lxy_partition        => false,
        :gmsh_filename        => "./problems/MHD/fluxEmergenceLeake2014/" * mesh,
        :lfilter              => false,
        :outformat            => "vtk",
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :output_dir           => "./output/",
        :loutput_pert         => false,
        :linitial_refine      => false,
        :ladapt               => false,
    ) #Dict
    haskey(ENV, "LK_RESTART_TIME") && (inputs[:lk_restart_time] = parse(Float64, ENV["LK_RESTART_TIME"]))

    return inputs
end
