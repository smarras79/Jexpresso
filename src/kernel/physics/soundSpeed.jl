function soundSpeed(npoin, mp, p_m, neqs, integrator, SD, ::TOTAL)
    
    # Jexpresso's own communicator, NOT MPI.COMM_WORLD. Under MPMD coupling
    # COMM_WORLD also contains Alya's ranks, which never call soundSpeed, so a
    # collective on it deadlocks every Jexpresso rank at the first diagnostic
    # output. get_mpi_comm() returns COMM_WORLD in standalone runs, so this is
    # a no-op there.
    comm = get_mpi_comm()
    # Physical constants
    PhysConst = PhysicalConst{Float32}()
    pos::TInt = 2
    if (SD == NSD_2D())
        pos = 3
    elseif (SD == NSD_3D())
        pos = 4
    end
    # Initialize arrays
    ρ = integrator.u[1:npoin]
    if (size(mp.Tabs,1)>1)
        pos_p = neqs+1
        Tabs = mp.Tabs[1:npoin]
        p = p_m
    elseif !ENERGY_EQUATION_THETA[]
        # Total-energy form: slot pos+1 holds ρE, so the ρθ gas law does
        # not apply; p = (γ-1)(ρE - ½ρ|u|²) summed over the momentum slots.
        ρE = integrator.u[pos*npoin+1:(pos+1)*npoin]
        ke = zero(ρ)
        for m = 2:pos
            ρum = @view integrator.u[(m-1)*npoin+1:m*npoin]
            ke .+= ρum.^2
        end
        p = PhysConst.γm1 .* (ρE .- 0.5 .* ke ./ ρ)
    else
        # Slot pos+1 holds ρθ, and the gas law p = C0·(ρθ)^γ takes ρθ
        # DIRECTLY. Passing (ρ, ρθ) to perfectGasLaw_ρθtoP, as this line used
        # to, evaluates C0·(ρ·ρθ)^γ instead.
        ρθ = integrator.u[pos*npoin+1:(pos+1)*npoin]
        p  = PhysConst.C0 .* max.(ρθ, zero(eltype(ρθ))).^PhysConst.γ
    end

    # Compute speed of sound using vectorized operation

    c = sqrt.(PhysConst.γ .* p ./ ρ)

    # Find the maximum speed of sound
    max_c = MPI.Allreduce(maximum(c), MPI.MAX, comm)

    return max_c
end

function soundSpeed(npoin, mp, p_m, neqs, integrator, SD, ::PERT)

    # Jexpresso's own communicator, NOT MPI.COMM_WORLD. Under MPMD coupling
    # COMM_WORLD also contains Alya's ranks, which never call soundSpeed, so a
    # collective on it deadlocks every Jexpresso rank at the first diagnostic
    # output. get_mpi_comm() returns COMM_WORLD in standalone runs, so this is
    # a no-op there.
    comm = get_mpi_comm()
    # Physical constants
    PhysConst = PhysicalConst{Float32}()
    pos::TInt = 2
    if (SD == NSD_2D())
        pos = 3
    elseif (SD == NSD_3D())
        pos = 4
    end
    # Initialize arrays
    ρ = integrator.u[1:npoin] + integrator.p.qp.qe[1:npoin]
    if (size(mp.Tabs,1)>1)
        pos_p = neqs+1
        Tabs = mp.Tabs[1:npoin]
        p = p_m
    else
        θ = (integrator.u[pos*npoin+1:(pos+1)*npoin] + integrator.p.qp.qe[pos*npoin+1:(pos+1)*npoin])./ρ
        # Compute pressure using vectorized operation
        p = perfectGasLaw_ρθtoP(PhysConst, ρ, θ)
    end

    # Compute speed of sound using vectorized operation

    c = sqrt.(PhysConst.γ .* p ./ ρ)

    # Find the maximum speed of sound
    max_c = MPI.Allreduce(maximum(c), MPI.MAX, comm)

    return max_c
end


function soundSpeed(npoin, mp, p_m, neqs, integrator, SD, ::THETA)
    
    # Physical constants
    PhysConst = PhysicalConst{Float32}()
    pos::TInt = 2
    if (SD == NSD_2D())
        pos = 3
    elseif (SD == NSD_3D())
        pos = 4
    end
    # Initialize arrays
    ρ = integrator.u[1:npoin]
    if (size(mp.Tabs,1)>1)
        pos_p = neqs+1
        Tabs = mp.Tabs[1:npoin]
        p = p_m
    else
        # p = C0·(ρθ)^γ takes ρθ directly — see the note in the TOTAL method.
        ρθ = integrator.u[pos*npoin+1:(pos+1)*npoin]
        p  = PhysConst.C0 .* max.(ρθ, zero(eltype(ρθ))).^PhysConst.γ
    end

    # Compute speed of sound using vectorized operation

    c = sqrt.(PhysConst.γ .* p ./ ρ)

    # Find the maximum speed of sound
    max_c = maximum(c)

    return max_c
end

# =============================================================================
#  CFL diagnostics
#
#  Each number below is a MAXIMUM over the mesh of a physical speed (or of a
#  diffusivity) times Δt over the smallest distance between two adjacent LGL
#  nodes.  Three things used to be wrong here, and every one of them biased
#  the printed CFL LOW — i.e. towards "there is plenty of margin left" on a
#  run that is in fact sitting near its stability limit:
#
#    1. the advective CFL was built from max(ρu, ρv), a MOMENTUM rather than
#       a velocity, and from the larger COMPONENT rather than from the speed
#       |u| = sqrt(u² + v²).  On CompEuler the error is a factor ρ, which is
#       ≈ 1.2 in the free stream of a sea-level case and 5-6 behind a strong
#       shock — so the number was not even wrong by a fixed factor;
#    2. the acoustic CFL was built from c alone.  The characteristic speed
#       of the Euler system is |u| + c; printing |u| and c in two separate
#       lines and taking the larger of the two understates the real number
#       by (|u|+c)/max(|u|,c), which is 4/3 at Mach 3 and 2 at Mach 1;
#    3. the length scale was Δelem/nop.  LGL nodes cluster towards the
#       element edges, so the smallest node gap is 0.69·Δelem/nop at nop = 4
#       and keeps shrinking like 1/nop² as the order goes up.  That gap is
#       what an explicit step has to resolve, and it is now what mesh.Δnode_s
#       carries (kernel/mesh/mesh.jl); Δelem/nop is still printed alongside
#       it so old logs remain comparable.
#
#  The viscous number had a fourth and larger problem: it was computed from
#  inputs[:μ], which for :visc_model => DSGS() is a dimensionless PER-EQUATION
#  MULTIPLIER and not a viscosity at all.  On CompEuler/ffs_step, with
#  :μ => [1.0, 4.0, 4.0, 4.0], the line printed 4.0·Δt/Δs² and had nothing to
#  do with the viscosity the run applies.  It now reads the DynSGS field that
#  the RHS actually assembled (params.μ_dsgs_pnode), converted to a KINEMATIC
#  diffusivity, because that is what sets the parabolic step limit.
# =============================================================================

# -----------------------------------------------------------------------------
# Node-wise maxima of |u|, of c and of |u| + c on THIS rank.
#
# Every system that reaches computeCFL stores CONSERVED variables: slot 1 is
# ρ (H for shallow water) and slots 2..nsd+1 are the momenta.  The speed is
# therefore the momentum divided by slot 1, node by node — never the maximum
# of the momentum field, and never the maxima of ρ and of ρu taken over
# different nodes and divided.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# PER-NODE LGL SPACING.
#
# The CFL numbers used to be built from two GLOBAL extrema: the largest wave
# speed ANYWHERE in the mesh over the smallest node spacing ANYWHERE in the
# mesh. On a uniform grid those live at the same place and the number is
# right. On a graded one they do not, and the print pairs a speed from the
# coarse far field with a length from the wall — on the 34x-graded
# shock_circle_M7 cylinder that turned a true parabolic number of 0.008 into a
# printed 0.312, a factor of 38, and the run was then "fixed" against a number
# that was never real.
#
# So each node gets the LGL spacing of the smallest element touching it, and
# the CFL maxima are taken over speed_i/Δ_i rather than max(speed)/min(Δ).
# The ratio Δnode_s/Δelem_s is the LGL fraction of the element side (0.17267 at
# nop = 4) and is the same for every element of a given order, so scaling each
# element's own size by it is exact on affine elements.
#
# Both global extrema are still PRINTED, because "max|u| = ..." and
# "max ν = ..." are useful on their own — they are simply no longer divided by
# a length from somewhere else.
# -----------------------------------------------------------------------------
function nodal_length_scale(mesh, SD)

    npoin = mesh.npoin
    Δ     = fill(Inf, npoin)

    Δe_s  = Float64(mesh.Δelem_s)
    Δn_s  = Float64(mesh.Δnode_s)
    ratio = (Δe_s > 0.0 && isfinite(Δn_s) && Δn_s > 0.0) ? Δn_s/Δe_s :
                                                           1.0/max(Float64(mesh.nop), 1.0)
    ngl = mesh.ngl
    @inbounds if SD == NSD_3D()
        for ie = 1:mesh.nelem
            Δe = Float64(mesh.Δelem[ie])*ratio
            (isfinite(Δe) && Δe > 0.0) || continue
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip = mesh.connijk[ie,i,j,k]
                (ip >= 1 && ip <= npoin) && (Δ[ip] = min(Δ[ip], Δe))
            end
        end
    else
        for ie = 1:mesh.nelem
            Δe = Float64(mesh.Δelem[ie])*ratio
            (isfinite(Δe) && Δe > 0.0) || continue
            for j = 1:ngl, i = 1:ngl
                ip = mesh.connijk[ie,i,j]
                (ip >= 1 && ip <= npoin) && (Δ[ip] = min(Δ[ip], Δe))
            end
        end
    end

    # Any node no element claimed (should not happen) falls back to the global
    # minimum, which is the old, conservative behaviour.
    fallback = (isfinite(Δn_s) && Δn_s > 0.0) ? Δn_s : 1.0
    @inbounds for ip = 1:npoin
        isfinite(Δ[ip]) || (Δ[ip] = fallback)
    end

    return Δ
end

function local_wave_speeds(npoin, neqs, mp, p_m, integrator, SD, Δnode)

    PhysConst = PhysicalConst{Float64}()

    nsd    = (SD == NSD_3D()) ? 3 : 2
    q      = integrator.u
    qe     = integrator.p.qp.qe
    lpert  = (integrator.p.SOL_VARS_TYPE == PERT())
    lmicro = size(mp.Tabs, 1) > 1
    lθ     = ENERGY_EQUATION_THETA[]

    # Slot carrying the energy variable (ρE or ρθ): 4 in 2D, 5 in 3D.  Systems
    # with fewer equations than that (shallow water) have no acoustic branch
    # of this kind and get c = 0.
    ien       = nsd + 2
    lacoustic = neqs >= ien
    tiny      = 1.0e-16
    # Shallow water (no acoustic slot): slot 1 is the depth H. The velocity
    # is Hu/H only down to the case's wet/dry threshold — the fluxes
    # desingularize it the same way — and the wave speed is √(gH). Dividing
    # by a thin film of 1e-16 printed |u| = 6e12 m/s and a CFL of 4e11 on a
    # run that was perfectly fine (SoliWaveIslandDSGS at t = 5 s, when the
    # wave reaches the island).
    lswe  = !lacoustic && neqs == nsd + 1
    g_swe = Float64(get(integrator.p.inputs, :dsgs_swe_g,    9.81))
    h_min = Float64(get(integrator.p.inputs, :dsgs_swe_hmin, 1.0e-3))

    velomax = 0.0
    cmax    = 0.0
    wavemax = 0.0
    # CFL per unit Δt, taken node by node against that node's OWN spacing.
    cflu_dt = 0.0
    cflc_dt = 0.0

    @inbounds for ip = 1:npoin

        ρ = lpert ? q[ip] + qe[ip] : q[ip]
        H = ρ
        ρ = max(ρ, tiny)   # a field that has already gone bad must still print

        ke = 0.0
        for m = 1:nsd
            ρum = lpert ? q[m*npoin + ip] + qe[m*npoin + ip] : q[m*npoin + ip]
            ke += ρum*ρum
        end
        vel = sqrt(ke)/ρ

        c = 0.0
        if lswe
            # the desingularized velocity of the shallow-water fluxes
            # (Kurganov & Petrova): Hu/H above the threshold, → 0 below it
            Hc  = max(H, 0.0)
            H4  = max(Hc, h_min)
            vel = sqrt(2.0)*Hc*sqrt(ke)/sqrt(Hc^4 + H4^4)
            c   = sqrt(g_swe*Hc)
        elseif lacoustic
            pl = 0.0
            if lmicro
                pl = p_m[ip]
            elseif !lθ
                # Total-energy form: p = (γ-1)(ρE - ½|ρu|²/ρ).
                ρE = lpert ? q[(ien-1)*npoin + ip] + qe[(ien-1)*npoin + ip] :
                             q[(ien-1)*npoin + ip]
                pl = PhysConst.γm1*(ρE - 0.5*ke/ρ)
            else
                # Slot ien holds ρθ, and the gas law is p = C0·(ρθ)^γ, so it
                # takes ρθ DIRECTLY.  The vectorised path this replaces called
                # perfectGasLaw_ρθtoP(PhysConst, ρ, ρθ) and so evaluated
                # C0·(ρ·ρθ)^γ — an extra factor ρ^γ inside the power, i.e. a
                # sound speed off by ρ^(γ/2) ≈ 1.14 at ρ = 1.2 and ≈ 0.75 at
                # the ρ = 0.4 of a stratosphere-topped column.
                ρθ = lpert ? q[(ien-1)*npoin + ip] + qe[(ien-1)*npoin + ip] :
                             q[(ien-1)*npoin + ip]
                pl = PhysConst.C0*(max(ρθ, 0.0))^PhysConst.γ
            end
            c = sqrt(max(PhysConst.γ*pl/ρ, 0.0))
        end

        velomax = max(velomax, vel)
        cmax    = max(cmax, c)
        wavemax = max(wavemax, vel + c)

        Δi      = Δnode[ip]
        cflu_dt = max(cflu_dt, vel/Δi)
        cflc_dt = max(cflc_dt, (vel + c)/Δi)
    end

    return velomax, cmax, wavemax, cflu_dt, cflc_dt
end

# -----------------------------------------------------------------------------
# Largest KINEMATIC diffusivity (m²/s) the viscous RHS actually applies on
# this rank.
#
# For DynSGS the viscosity is a field, recomputed every RHS call from the
# residual, so inputs[:μ] says nothing about its magnitude — it is only the
# per-equation multiplier that compute_dsgs_viscosity! has already folded into
# μ_dsgs.  Slot 1 of μ_dsgs_pnode is β = μ/‖ρ‖_{∞,K} and is already kinematic;
# every other slot is the DYNAMIC coefficient that multiplies a primitive
# gradient in a conserved-variable equation, so μ/ρ is the diffusivity that
# sets the parabolic limit.
# -----------------------------------------------------------------------------
function local_max_diffusivity(npoin, params, visc, Δnode)

    ldsgs = (params.VT == DSGS() || params.VT == DSGS_MHD() || params.VT == DSGS_SW()) &&
            size(params.μ_dsgs_pnode, 1) == npoin

    if !ldsgs
        # Constant-coefficient models (AV, Smagorinsky-with-fixed-μ, …): here
        # inputs[:μ] IS the coefficient, so the old behaviour is the right one.
        # ν is uniform, so the tightest node is simply the smallest one.
        νc = maximum(visc)
        Δm = minimum(Δnode)
        return νc, νc/(Δm*Δm)
    end

    q     = params.uaux
    qe    = params.qp.qe
    lpert = (params.SOL_VARS_TYPE == PERT())
    neqsν = size(params.μ_dsgs_pnode, 2)
    tiny  = 1.0e-16

    # Which slots hold a DYNAMIC coefficient (to be divided by ρ) and which a
    # KINEMATIC one (a diffusivity already). The MHD kernel stores the
    # magnetic and ψ slots (6-9) as kinematic resistivities, and with
    # :dsgs_nodal_rho or :dsgs_conserved its momentum/energy slots too (the
    # ρ factor is then applied per quadrature point at assembly, or absent).
    # Dividing those by the 7e-9 of a solar corona printed a "max ν" of 1e7
    # for a run whose real parabolic number was 0.04.
    mhd    = (params.VT == DSGS_MHD())
    # Euler kernels: the passive-tracer slots carry the KINEMATIC ν. They
    # start after the energy slot, which is 4 in 2D (ρ, ρu, ρv, ρθ) and 5 in
    # 3D (ρ, ρu, ρv, ρw, ρθ).
    euler   = (params.VT == DSGS())
    tracer0 = (params.SD == NSD_3D()) ? 6 : 5
    allkin = (mhd && (get(params.inputs, :dsgs_nodal_rho, false) || get(params.inputs, :dsgs_conserved, false))) ||
             params.VT == DSGS_SW()      # shallow water: one kinematic ν on (H, Hu, Hv)

    ν       = 0.0
    parab   = 0.0        # max over nodes of ν_i/Δ_i², i.e. the parabolic
                         # number per unit Δt
    @inbounds for ip = 1:npoin
        ρ = lpert ? q[ip,1] + qe[ip,1] : q[ip,1]
        ρ = max(ρ, tiny)
        νi = params.μ_dsgs_pnode[ip,1]                  # β / mass diffusion, kinematic
        for ieq = 2:neqsν
            kin = allkin || (mhd && ieq >= 6) || (euler && ieq >= tracer0)
            νi = max(νi, kin ? params.μ_dsgs_pnode[ip,ieq] : params.μ_dsgs_pnode[ip,ieq]/ρ)
        end
        ν     = max(ν, νi)
        Δi    = Δnode[ip]
        parab = max(parab, νi/(Δi*Δi))
    end

    return ν, parab
end

# -----------------------------------------------------------------------------
# One-time parabolic-number check for the DynSGS models, called after the
# warm-up step (TimeIntegrators.jl) once μ_dsgs_pnode holds the coefficient
# of the first step. A residual viscosity sits at its first-order cap
# C_max·Δ·(|u|+c) wherever the initial condition has a kink (the cone edge of
# a θ bubble, a tracer top-hat, a diaphragm), and an explicit RK step can only
# carry ν·Δt/Δx_min² up to O(1): CompEuler/thetaTracers with DSGS() blew up
# at the first step at 1.6 while its deck's SMAG run reads 2e-4 for the same
# Δt. The diagnostics callback prints the same number, but only at the first
# output time, which the run never reached. Nothing is changed here; the
# warning names the number and the two knobs (Δt, the :μ multipliers).
# -----------------------------------------------------------------------------
function dsgs_first_step_check(params, inputs, SD::Union{NSD_2D, NSD_3D})
    ldsgs = (params.VT == DSGS() || params.VT == DSGS_MHD() || params.VT == DSGS_SW()) &&
            size(params.μ_dsgs_pnode, 1) == params.mesh.npoin
    ldsgs || return nothing
    comm    = get_mpi_comm()
    # Per-node, like computeCFL: max_i(ν_i/Δ_i²), not max(ν)/min(Δ)². On a
    # graded mesh the old pairing took ν from the coarse far field and Δ from
    # the wall and overstated the number by the square of the grading.
    Δnode_v = nodal_length_scale(params.mesh, SD)
    νmax_l, parab_dt_l = local_max_diffusivity(params.mesh.npoin, params, inputs[:μ], Δnode_v)
    buf     = MPI.Allreduce([νmax_l, parab_dt_l], MPI.MAX, comm)
    νmax, parab_dt = buf[1], buf[2]
    pnum    = parab_dt*Float64(inputs[:Δt])
    if pnum > 0.5 && MPI.Comm_rank(comm) == 0
        @warn @sprintf("DynSGS after the first step: max ν = %.3e m²/s, max(ν·Δt/Δx²) over the nodes = %.2f (Δt = %g s). Above ~0.5 the explicit step cannot carry the diffusion and the run blows up: reduce :Δt (or the :μ multipliers of the slots that carry the largest coefficient).",
                       νmax, pnum, Float64(inputs[:Δt]))
    end
    return nothing
end
dsgs_first_step_check(params, inputs, SD) = nothing

function computeCFL(npoin, neqs, mp, p, dt, Δs, integrator, SD::NSD_1D; visc=[0.0])
    nothing
end

function computeCFL(npoin, neqs, mp, p, dt, Δs, integrator, SD::Union{NSD_2D, NSD_3D}; visc=[0.0])

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    nsd = (SD == NSD_3D()) ? 3 : 2
    if size(integrator.u, 1) < (nsd+1)*npoin
        return nothing
    end

    Δnode_v = nodal_length_scale(integrator.p.mesh, SD)

    velomax_l, cmax_l, wavemax_l, cflu_dt_l, cflc_dt_l =
        local_wave_speeds(npoin, neqs, mp, p, integrator, SD, Δnode_v)
    νmax_l, parab_dt_l = local_max_diffusivity(npoin, integrator.p, visc, Δnode_v)

    # One packed Allreduce instead of four latency-bound round trips.  Doing
    # the reduction here rather than inside soundSpeed() is also what makes
    # the printed numbers partition-independent: a rank-local maximum would
    # make the diagnostic depend on how the mesh happens to be split.
    local_buf  = [velomax_l, cmax_l, wavemax_l, νmax_l,
                  cflu_dt_l, cflc_dt_l, parab_dt_l]
    global_buf = MPI.Allreduce(local_buf, MPI.MAX, comm)
    velomax, cmax, wavemax, νmax =
        global_buf[1], global_buf[2], global_buf[3], global_buf[4]
    cflu_dt, cflc_dt, parab_dt =
        global_buf[5], global_buf[6], global_buf[7]

    # Smallest gap between two adjacent LGL nodes anywhere in the mesh; it
    # already reflects AMR, because it is measured on the refined grid.
    # Fall back to the Δelem/nop that the caller passes if a mesh path has
    # not filled it (Δnode_s == 0).
    Δnode = Float64(integrator.p.mesh.Δnode_s)
    Δ     = (isfinite(Δnode) && Δnode > 0.0) ? Δnode : Float64(Δs)
    dtf   = Float64(dt)

    # Each of these is now max_i(speed_i/Δ_i)·Δt, NOT max(speed)·Δt/min(Δ):
    # the speed and the length come from the same node. See the comment on
    # nodal_length_scale above for what the old pairing cost.
    cfl_u    = cflu_dt*dtf             # advective
    cfl_c    = cflc_dt*dtf             # acoustic: the |u| + c characteristic
    cfl_visc = parab_dt*dtf            # parabolic

    println_rank(@sprintf(" #  Δx_min (LGL) : %.4e m   (Δelem/nop = %.4e m)",
                          Δ, Float64(Δs)); msg_rank = rank)
    println_rank(@sprintf(" #  Advective CFL: %.6f   max|u|      = %.4e m/s",
                          cfl_u, velomax); msg_rank = rank)
    println_rank(@sprintf(" #  Acoustic  CFL: %.6f   max(|u|+c)  = %.4e m/s  (max c = %.4e)",
                          cfl_c, wavemax, cmax); msg_rank = rank)
    println_rank(@sprintf(" #  Viscous   CFL: %.6f   max ν       = %.4e m²/s",
                          cfl_visc, νmax); msg_rank = rank)

    return nothing
end
