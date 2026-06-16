function soundSpeed(npoin, mp, p_m, neqs, integrator, SD, ::TOTAL)
    
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
        θ = integrator.u[pos*npoin+1:(pos+1)*npoin]
        # Compute pressure using vectorized operation
        p = perfectGasLaw_ρθtoP(PhysConst, ρ, θ)
    end

    # Compute speed of sound using vectorized operation

    c = sqrt.(PhysConst.γ .* p ./ ρ)

    # Find the maximum speed of sound
    max_c = maximum(c)

    return max_c
end

function soundSpeed(npoin, mp, p_m, neqs, integrator, SD, ::PERT)

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
    max_c = maximum(c)

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
	θ = (integrator.u[pos*npoin+1:(pos+1)*npoin])
        # Compute pressure using vectorized operation
        p = perfectGasLaw_ρθtoP(PhysConst, ρ, θ)
    end

    # Compute speed of sound using vectorized operation

    c = sqrt.(PhysConst.γ .* p ./ ρ)

    # Find the maximum speed of sound
    max_c = maximum(c)

    return max_c
end

function computeCFL(npoin, neqs, mp, p, dt, Δs, integrator, SD::NSD_1D; visc=[0.0])
    nothing
end

function computeCFL(npoin, neqs, mp, p, dt, Δs, integrator, SD::NSD_2D; visc=[0.0])

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    if size(integrator.u)[1] >= 3*npoin
        # Compute local maxima first, then do ONE Allreduce on a packed
        # vector instead of 2-3 sequential MAX reductions. Each Allreduce
        # is a latency-bound round trip; batching collapses them.
        # Slots: [umax, vmax, cmax]. cmax local = 0 when neqs<4 -> ignored.
        if (integrator.p.SOL_VARS_TYPE == PERT())
            idx  = npoin
            umax_local = maximum(integrator.u[idx+1:2*npoin] + integrator.p.qp.qe[idx+1:2*npoin])
            idx  = 2*npoin
            vmax_local = maximum(integrator.u[idx+1:3*npoin] + integrator.p.qp.qe[idx+1:3*npoin])
        else
            idx  = npoin
            umax_local = maximum(integrator.u[idx+1:2*npoin])
            idx  = 2*npoin
            vmax_local = maximum(integrator.u[idx+1:3*npoin])
        end

        # speed of sound (only for systems with >= 4 equations, e.g. CompEuler)
        # Note: previously soundSpeed returned a per-rank-local max and the
        # printed acoustic CFL was wrong on multi-rank runs. Folding `cmax`
        # into the same Allreduce fixes that as a side effect.
        cmax_local = (neqs >= 4) ?
            soundSpeed(npoin, mp, p, neqs, integrator, SD, integrator.p.SOL_VARS_TYPE) :
            zero(eltype(integrator.u))

        local_buf  = [umax_local, vmax_local, cmax_local]
        global_buf = MPI.Allreduce(local_buf, MPI.MAX, comm)
        umax, vmax, cmax = global_buf[1], global_buf[2], global_buf[3]

        velomax = max(umax, vmax)
        cfl_u   = velomax*dt/Δs        # Advective CFL
        cfl_c   = (neqs >= 4) ? cmax*dt/Δs : 0.0   # Acoustic CFL

        Δs2      = Δs*Δs
        μ        = maximum(visc)
        λ        = 1.0 #free parameter
        cfl_visc = dt*λ*μ/Δs2 #Viscous CFL

        println_rank(" #  Advective CFL: ", cfl_u;    msg_rank = rank) #, suppress = mesh.msg_suppress)
        println_rank(" #  Acoustic  CFL: ", cfl_c;    msg_rank = rank) #, suppress = mesh.msg_suppress)
        println_rank(" #  Viscous   CFL: ", cfl_visc; msg_rank = rank) #, suppress = mesh.msg_suppress)
    else
        nothing
    end

end

function computeCFL(npoin, neqs, mp, p, dt, Δs, integrator, SD::NSD_3D; visc=[0.0])

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    if size(integrator.u)[1] >= 4*npoin
        # Compute local maxima first, then do ONE Allreduce on a packed
        # vector instead of 4 sequential MAX reductions. Each Allreduce
        # is a latency-bound round trip; batching collapses them.
        # Slots: [umax, vmax, wmax, cmax, μmax].
        if (integrator.p.SOL_VARS_TYPE == PERT())
            idx  = npoin
            umax_local = maximum(integrator.u[idx+1:2*npoin] + integrator.p.qp.qe[idx+1:2*npoin])
            idx  = 2*npoin
            vmax_local = maximum(integrator.u[idx+1:3*npoin] + integrator.p.qp.qe[idx+1:3*npoin])
            idx  = 3*npoin
            wmax_local = maximum(integrator.u[idx+1:4*npoin] + integrator.p.qp.qe[idx+1:4*npoin])
        else
            idx  = npoin
            umax_local = maximum(integrator.u[idx+1:2*npoin])
            idx  = 2*npoin
            vmax_local = maximum(integrator.u[idx+1:3*npoin])
            idx  = 3*npoin
            wmax_local = maximum(integrator.u[idx+1:4*npoin])
        end

        # speed of sound. soundSpeed previously returned a per-rank-local
        # max and the printed acoustic CFL was wrong on multi-rank runs;
        # folding cmax into the same Allreduce fixes that as a side effect.
        cmax_local = soundSpeed(npoin, mp, p, neqs, integrator, SD, integrator.p.SOL_VARS_TYPE)
        μmax_local = maximum(integrator.p.μ_max)

        local_buf  = [umax_local, vmax_local, wmax_local, cmax_local, μmax_local]
        global_buf = MPI.Allreduce(local_buf, MPI.MAX, comm)
        umax, vmax, wmax, cmax, μmax =
            global_buf[1], global_buf[2], global_buf[3], global_buf[4], global_buf[5]

        velomax = max(umax, vmax, wmax)
        cfl_u   = velomax*dt/Δs   # Advective CFL
        cfl_c   = cmax*dt/Δs      # Acoustic CFL

        Δs2      = Δs*Δs
        μ        = max(maximum(visc), μmax)
        λ        = 2.0 #free parameter
        cfl_visc = dt*λ*μ/Δs2 #Viscous CFL

        
        println_rank(" #  Advective CFL: ", cfl_u;    msg_rank = rank)
        println_rank(" #  Acoustic  CFL: ", cfl_c;    msg_rank = rank)
        println_rank(" #  Viscous   CFL: ", cfl_visc; msg_rank = rank) #, suppress = mesh.msg_suppress)
    end
end
