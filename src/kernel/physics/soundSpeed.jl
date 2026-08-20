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
function local_wave_speeds(npoin, neqs, mp, p_m, integrator, SD)

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

    velomax = 0.0
    cmax    = 0.0
    wavemax = 0.0

    @inbounds for ip = 1:npoin

        ρ = lpert ? q[ip] + qe[ip] : q[ip]
        ρ = max(ρ, tiny)   # a field that has already gone bad must still print

        ke = 0.0
        for m = 1:nsd
            ρum = lpert ? q[m*npoin + ip] + qe[m*npoin + ip] : q[m*npoin + ip]
            ke += ρum*ρum
        end
        vel = sqrt(ke)/ρ

        c = 0.0
        if lacoustic
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
    end

    return velomax, cmax, wavemax
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
function local_max_diffusivity(npoin, params, visc)

    ldsgs = (params.VT == DSGS() || params.VT == DSGS_MHD()) &&
            size(params.μ_dsgs_pnode, 1) == npoin

    if !ldsgs
        # Constant-coefficient models (AV, Smagorinsky-with-fixed-μ, …): here
        # inputs[:μ] IS the coefficient, so the old behaviour is the right one.
        return maximum(visc)
    end

    q     = params.uaux
    qe    = params.qp.qe
    lpert = (params.SOL_VARS_TYPE == PERT())
    neqsν = size(params.μ_dsgs_pnode, 2)
    tiny  = 1.0e-16

    ν = 0.0
    @inbounds for ip = 1:npoin
        ρ = lpert ? q[ip,1] + qe[ip,1] : q[ip,1]
        ρ = max(ρ, tiny)
        ν = max(ν, params.μ_dsgs_pnode[ip,1])           # β, already kinematic
        for ieq = 2:neqsν
            ν = max(ν, params.μ_dsgs_pnode[ip,ieq]/ρ)   # dynamic → kinematic
        end
    end

    return ν
end

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

    velomax_l, cmax_l, wavemax_l = local_wave_speeds(npoin, neqs, mp, p, integrator, SD)
    νmax_l                       = local_max_diffusivity(npoin, integrator.p, visc)

    # One packed Allreduce instead of four latency-bound round trips.  Doing
    # the reduction here rather than inside soundSpeed() is also what makes
    # the printed numbers partition-independent: a rank-local maximum would
    # make the diagnostic depend on how the mesh happens to be split.
    local_buf  = [velomax_l, cmax_l, wavemax_l, νmax_l]
    global_buf = MPI.Allreduce(local_buf, MPI.MAX, comm)
    velomax, cmax, wavemax, νmax =
        global_buf[1], global_buf[2], global_buf[3], global_buf[4]

    # Smallest gap between two adjacent LGL nodes anywhere in the mesh; it
    # already reflects AMR, because it is measured on the refined grid.
    # Fall back to the Δelem/nop that the caller passes if a mesh path has
    # not filled it (Δnode_s == 0).
    Δnode = Float64(integrator.p.mesh.Δnode_s)
    Δ     = (isfinite(Δnode) && Δnode > 0.0) ? Δnode : Float64(Δs)
    dtf   = Float64(dt)

    cfl_u    = velomax*dtf/Δ           # advective
    cfl_c    = wavemax*dtf/Δ           # acoustic: the |u| + c characteristic
    cfl_visc = νmax*dtf/(Δ*Δ)          # parabolic

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
