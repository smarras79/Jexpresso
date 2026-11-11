#---------------------------------------------------------------------------------
# 2.5D analogue of the flux-emergence study of Leake, Linton & Antiochos, "The role
# of reconnection in the onset of solar eruptions", ApJ 934:10 (2022): the plane
# perpendicular to the rope, invariant along it. See README.md.
#
#   julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022_2D
#   mpiexec -n 14 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022_2D
#
# Environment overrides (README): LK2D_REF (mesh), LK2D_TEND, LK2D_OUT_DT, LK2D_DIAG,
# LK2D_RESTART_DIR (+ LK2D_RESTART_TIME).
#   LK2D_REF 1:  32 x 48 elements   (the 3D LK22_REF = 2 cross-section; default)
#            2:  64 x 96 elements   (1/2 of the paper's resolution at :nop => 4)
#            3: 128 x 192 elements  (the paper's 512 x 768 cells at :nop => 4)
#---------------------------------------------------------------------------------

# Positivity floors as the Runge-Kutta stage limiter (as the 3D Leake cases): ρ ≥ f ρ_e(h)
# and p ≥ f p_e(h) of the reference atmosphere, p raised through E at fixed velocity and
# field. Not conservative where it fires; the hits are counted over all ranks and reported
# by rank 0 each time the total passes a power of ten.
if !@isdefined(lk2d_floor_frac)
    const lk2d_floor_frac = Ref{Float64}(0.2)
    const lk2d_floor_hits = Ref{Int}(0)
    const lk2d_floor_next = Ref{Int}(1)
end

# CFL monitor (the time step is fixed, the wave speeds are not): every 100 stages the largest
# (|v| + c_f) Δt/Δh over the nodes, with c_f = √((γp + |B|²)/ρ) and Δh the smallest LGL gap of
# the elements around the node, is reported by rank 0 when it has grown by 20 %, as a WARNING
# above 0.5 (the step is set for ≈ 0.25 at |v| + c_f = 80: README).
if !@isdefined(lk2d_cfl_dh)
    const lk2d_cfl_dh    = Ref{Vector{Float64}}(Float64[])
    const lk2d_cfl_count = Ref{Int}(0)
    const lk2d_cfl_shown = Ref{Float64}(0.0)
    const lk2d_diag_calls = Ref{Int}(0)
end

function lk2d_node_spacing(mesh)
    ngl = Int(mesh.ngl)
    h   = fill(Inf, Int(mesh.npoin))
    for ie = 1:Int(mesh.nelem)
        he = Inf
        for m = 1:ngl-1, (a, b) in ((mesh.connijk[ie,m,1], mesh.connijk[ie,m+1,1]),
                                    (mesh.connijk[ie,1,m], mesh.connijk[ie,1,m+1]))
            he = min(he, sqrt((mesh.coords[1,b] - mesh.coords[1,a])^2 + (mesh.coords[2,b] - mesh.coords[2,a])^2))
        end
        for j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[ie,i,j]
            h[ip] = min(h[ip], he)
        end
    end
    return h
end

function lk2d_positivity_limiter!(u, integrator, p, t)
    npoin = p.mesh.npoin
    qe    = p.qp.qe
    frac  = lk2d_floor_frac[]
    γm1   = γ_mhd - 1.0
    nhit  = 0
    lk2d_cfl_count[] += 1
    docfl = lk2d_cfl_count[] % 100 == 1
    if docfl && length(lk2d_cfl_dh[]) != npoin
        lk2d_cfl_dh[] = lk2d_node_spacing(p.mesh)
    end
    dh  = lk2d_cfl_dh[]
    Δt  = Float64(p.Δt)
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
    comm = MPI.COMM_WORLD
    if docfl
        g = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(cfl, MPI.MAX, comm) : cfl
        if g > 1.2*lk2d_cfl_shown[]
            if MPI.Comm_rank(comm) == 0       # the |v|, c_f printed are those of rank 0's largest
                print(string(g > 0.5 ? " # WARNING CFL: " : " # CFL: ", "max (|v| + c_f) Δt/Δh = ", round(g, sigdigits = 3),
                             " at t = ", round(t, digits = 3), " (rank 0: |v| = ", round(vmx, sigdigits = 3), ", c_f = ",
                             round(cmx, sigdigits = 3), ")", g > 0.5 ? "; reduce :Δt" : "", "\n"))
                flush(stdout)
            end
            lk2d_cfl_shown[] = g
        end
    end
    # One report for all ranks, from rank 0, when the total passes 1, 10, 100, ...: every rank
    # calls the stage limiter at every Runge-Kutta stage, so the sum is a safe collective. t is
    # the stage time t_n + c_i Δt (in t₀), not an output time.
    nall = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(nhit, +, comm) : nhit
    lk2d_floor_hits[] += nall
    if nall > 0 && lk2d_floor_hits[] >= lk2d_floor_next[]
        if MPI.Comm_rank(comm) == 0
            print(string(" # positivity floors: ", lk2d_floor_hits[], " node hits so far on all ranks (",
                         nall, " in this stage) at t = ", round(t, digits = 3), "\n"))
            flush(stdout)
        end
        while lk2d_floor_next[] <= lk2d_floor_hits[]
            lk2d_floor_next[] *= 10
        end
    end
    lk2d_diag!(u, p, t)
    return nothing
end

# Stage diagnostics (LK2D_DIAG=n in the environment: every n-th stage; off by default), as the
# 2014 case's LK_DIAG: the largest |v| and |v| + c_f, the smallest ρ/ρ_e, the largest and
# smallest p/p_e, T, the nodes at the floors and the heights they span, the largest DynSGS +
# explicit coefficient (ρ, momentum, B slots) and the explicit-diffusion number μΔt/Δh², then
# where each extreme is.
function lk2d_diag!(u, p, t)
    n = parse(Int, get(ENV, "LK2D_DIAG", "0"))
    n <= 0 && return nothing
    lk2d_diag_calls[] += 1
    lk2d_diag_calls[] % n == 0 || return nothing
    npoin = p.mesh.npoin; qe = p.qp.qe; γm1 = γ_mhd - 1.0; X = p.mesh.coords
    vmax = 0.0; wmax = 0.0; ρrmin = Inf; prmax = 0.0; prmin = Inf; Tmax = 0.0
    nρfl = 0; npfl = 0; hlo = Inf; hhi = -Inf
    xw = (0.0, 0.0); xv = xw; xρ = xw; xp = xw; xpm = xw
    frac = 1.000001*lk2d_floor_frac[]
    @inbounds for ip = 1:npoin
        ρ = u[ip]; ru = u[npoin+ip]; rv = u[2npoin+ip]; E = u[3npoin+ip]; rw = u[4npoin+ip]
        Bx = u[5npoin+ip]; By = u[6npoin+ip]; Bz = u[7npoin+ip]; ψ = u[8npoin+ip]
        B2 = Bx*Bx + By*By + Bz*Bz
        v  = sqrt(ru*ru + rv*rv + rw*rw)/ρ
        pr = γm1*(E - 0.5*(ru*ru + rv*rv + rw*rw)/ρ - 0.5*B2 - 0.5*ψ*ψ)
        cf = sqrt(max((γ_mhd*pr + B2)/ρ, 0.0))
        x  = (X[1,ip], X[2,ip])
        v > vmax && (vmax = v; xv = x)
        v + cf > wmax && (wmax = v + cf; xw = x)
        rr = ρ/qe[ip,1]; pp = pr/qe[ip,end]
        rr < ρrmin && (ρrmin = rr; xρ = x)
        pp > prmax && (prmax = pp; xp = x)
        pp < prmin && (prmin = pp; xpm = x)
        rr <= frac && (nρfl += 1)
        pp <= frac && (npfl += 1)
        (rr <= frac || pp <= frac) && (hlo = min(hlo, x[2]); hhi = max(hhi, x[2]))
        Tmax = max(Tmax, lk2d_mf[]*pr/ρ)
    end
    μd = p.μ_dsgs; Δt = Float64(p.Δt); dh = lk2d_cfl_dh[]
    μρ = 0.0; μm = 0.0; μB = 0.0; Dmax = 0.0; xD = (0.0, 0.0)
    if size(μd, 1) == p.mesh.nelem && length(dh) == npoin
        @inbounds for ie = 1:p.mesh.nelem
            μρ = max(μρ, μd[ie,1]); μm = max(μm, μd[ie,2]); μB = max(μB, μd[ie,6])
            ip = p.mesh.connijk[ie,2,2]
            D  = maximum(@view μd[ie,1:8])*Δt/dh[ip]^2
            D > Dmax && (Dmax = D; xD = (X[1,ip], X[2,ip]))
        end
    end
    comm = MPI.COMM_WORLD
    glo = MPI.Allreduce([vmax, wmax, -ρrmin, prmax, -prmin, Tmax, μρ, μm, μB, Dmax, -hlo, hhi], MPI.MAX, comm)
    cnt = MPI.Allreduce([nρfl, npfl], +, comm)
    if MPI.Comm_rank(comm) == 0
        println(" #DIAG t=", round(t, digits=4), " |v|max=", round(glo[1], sigdigits=4),
                " (|v|+cf)max=", round(glo[2], sigdigits=4), " ρ/ρe min=", round(-glo[3], sigdigits=3),
                " p/pe max=", round(glo[4], sigdigits=4), " p/pe min=", round(-glo[5], sigdigits=3),
                " Tmax=", round(glo[6], sigdigits=4), " nodes at ρ/p floor=", cnt[1], "/", cnt[2],
                " μmax(ρ,m,B)=", round(glo[7], sigdigits=3), ",", round(glo[8], sigdigits=3), ",",
                round(glo[9], sigdigits=3), " μΔt/Δh²max=", round(glo[10], sigdigits=3),
                cnt[1] + cnt[2] > 0 ? string(" floors at y=", round(-glo[11], digits=1), "..", round(glo[12], digits=1)) : "")
    end
    for (val, g, x, name) in ((wmax, glo[2], xw, "(|v|+cf)max"), (vmax, glo[1], xv, "|v|max"),
                              (-ρrmin, glo[3], xρ, "ρ/ρe min"), (prmax, glo[4], xp, "p/pe max"),
                              (-prmin, glo[5], xpm, "p/pe min"), (Dmax, glo[10], xD, "μΔt/Δh² max"))
        val == g && println(" #DIAG   ", name, " at (x, y) = ", round.(x, digits=2), " (rank ", MPI.Comm_rank(comm), ")")
    end
    flush(stdout)
    return nothing
end

# Explicit dissipation of the paper (Sec. 2.1: "both equal to 0.01" in normalized units,
# Re = Rm = 100): kinematic ν on the momentum and uniform η on the field, per element
# centroid, added to the DynSGS coefficient (physics, not stabilization: README). In 2D the
# function takes (x, y).
lk2d_dissipation(x, y) = (0.01, 0.01)

function user_inputs()

    LK2D_REF = parse(Int, get(ENV, "LK2D_REF", "1"))
    mesh     = ("LK22_2D_32x48.msh", "LK22_2D_64x96.msh", "LK22_2D_128x192.msh")[LK2D_REF]
    # Δt for (|v| + c_f) Δt/Δh ≈ 0.25 at |v| + c_f ≈ 80 on the smallest LGL gap Δh of each mesh
    # (0.173 of the 4, 2, 1 L₀ tall bottom elements), as the 3D deck at the same element size.
    Δt       = (2.0e-3, 1.0e-3, 5.0e-4)[LK2D_REF]
    tend     = parse(Float64, get(ENV, "LK2D_TEND", "400.0"))
    dtout    = parse(Float64, get(ENV, "LK2D_OUT_DT", "5.0"))

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(stage_limiter! = lk2d_positivity_limiter!,
                                                      williamson_condition = false),
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => tend,       # 400: the 3D deck's (its simulation 4 erupts at 150-200 t₀)
        :diagnostics_at_times => filter(<=(tend), sort(unique(vcat(collect(0.0:dtout:tend), [142.0, 184.0])))),
        :conservation_every   => 10,
        :conservation_slots   => [1, 4],
        :lstep_heartbeat      => true,       # "#   step N   t = …" every 100 steps
        :restart_time         => 0.0,
        :lrestart             => false,
        :lk2d_restart_dir     => get(ENV, "LK2D_RESTART_DIR", ""),   # "" = start from the initial state
        :lsource              => true,       # gravity, ψ damping, velocity damping
        :lwell_balanced       => true,       # subtract the discrete RHS of qe (atmosphere + dipole)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # The problem (README):
        #---------------------------------------------------------------------------
        :lk2d_orientation     => "antiparallel", # the 3D θ = π (simulation 4, erupts in 3D); "parallel": θ = 0 (stable)
        :lk2d_dipole          => true,       # false: no coronal field
        :lk2d_rope            => true,       # false: the background alone, which must stay at rest
        :lk2d_D               => 150.0,      # line-dipole strength: 19.5 G at (0, 0), the 3D dipole's value
        :lk2d_hd              => -100.0,     # line-dipole depth
        :lk2d_buoyancy        => 1.0,        # b: 1 = the middle of the 3D rope (isothermal), 0 = neutrally buoyant
        :lk2d_cz_gradient     => 1.0,        # Eq. 23 as printed (exactly adiabatic with m_f = 1)
        :lk2d_mf              => 1.0,        # p = ρT/m_f (README of the 3D case)
        :lk2d_glm_alpha       => 1.0,        # ψ damping rate α_p c_h/Δh
        :lk2d_damp_x_frac     => 0.75,       # side damping from |x| = 3L_x/8 (Eq. 21)
        :lk2d_damp_z_frac     => 0.875,      # top damping from y = 7 max(y)/8 (Eq. 22)
        :mhd_dissipation      => lk2d_dissipation,   # the paper's ν = η = 0.01; remove for ideal MHD + DynSGS
        #---------------------------------------------------------------------------
        # DynSGS in conserved form on the departure from the reference state, with the
        # constants of the 3D Leake decks.
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
        :dsgs_fast_floors     => true,
        :dsgs_thermal_residual => true,
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh: LK22_2D.geo (the paper's stretching in the plane, element boundaries every
        # 16, 8, 4 cells)
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/fluxEmergenceLeake2022_2D/" * mesh,
        # Balanced (METIS) partition: the default column partition bins the bounding box
        # uniformly, which leaves ranks empty on this stretched mesh (3D case's README).
        :lxy_partition        => false,
        :lfilter              => false,
        #---------------------------------------------------------------------------
        # Output: VTK (ParaView, tools/plot_fluxemergence_vtu.py, the restart), with the
        # extra fields log10_ρ, log10_p, log10_β.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :plot_log10           => ["ρ", "p", "β"],
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :output_dir           => "./output/",
        :loutput_pert         => false,
        :linitial_refine      => false,
        :ladapt               => false,
    ) #Dict
    haskey(ENV, "LK2D_RESTART_TIME") && (inputs[:lk2d_restart_time] = parse(Float64, ENV["LK2D_RESTART_TIME"]))

    return inputs
end
