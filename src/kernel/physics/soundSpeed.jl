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
        θ = (integrator.u[pos*npoin+1:(pos+1)*npoin] + integrator.p.q.qe[pos*npoin+1:(pos+1)*npoin])./ρ
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

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    if size(integrator.u)[1] >= 3*npoin
        #u
        if (integrator.p.SOL_VARS_TYPE == PERT())
            idx  = npoin
            umax = MPI.Allreduce(maximum(integrator.u[idx+1:2*npoin] + integrator.p.qp.qe[idx+1:2*npoin]), MPI.MAX, comm)

            idx  = 2*npoin
            vmax = MPI.Allreduce(maximum(integrator.u[idx+1:3*npoin] + integrator.p.qp.qe[idx+1:3*npoin]), MPI.MAX, comm)
        else
            idx  = npoin
            umax = MPI.Allreduce(maximum(integrator.u[idx+1:2*npoin]), MPI.MAX, comm)
        
            idx  = 2*npoin
            vmax = MPI.Allreduce(maximum(integrator.u[idx+1:3*npoin]), MPI.MAX, comm)
        end
        velomax = max(umax, vmax)
        
        #speed of sound
        c     = soundSpeed(npoin, mp, p, neqs, integrator, SD, integrator.p.SOL_VARS_TYPE)
        cfl_u = velomax*dt/Δs #Advective CFL
        cfl_c = c*dt/Δs       #Acoustic CFL
        
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

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    if size(integrator.u)[1] >= 4*npoin
        #u
        if (integrator.p.SOL_VARS_TYPE == PERT())
            idx  = npoin
            umax = MPI.Allreduce(maximum(integrator.u[idx+1:2*npoin] + integrator.p.qp.qe[idx+1:2*npoin]), MPI.MAX, comm)

            idx  = 2*npoin
            vmax = MPI.Allreduce(maximum(integrator.u[idx+1:3*npoin] + integrator.p.qp.qe[idx+1:3*npoin]), MPI.MAX, comm)
            
            idx  = 3*npoin
            wmax = MPI.Allreduce(maximum(integrator.u[idx+1:4*npoin] + integrator.p.qp.qe[idx+1:4*npoin]), MPI.MAX, comm)
        else
            idx  = npoin
            umax = MPI.Allreduce(maximum(integrator.u[idx+1:2*npoin]), MPI.MAX, comm)

            idx  = 2*npoin
            vmax = MPI.Allreduce(maximum(integrator.u[idx+1:3*npoin]), MPI.MAX, comm)
            
            idx  = 3*npoin
            wmax = MPI.Allreduce(maximum(integrator.u[idx+1:4*npoin]), MPI.MAX, comm)
        end
        
        velomax = max(umax, vmax, wmax)
        
        #speed of sound
        c     = soundSpeed(npoin, mp, p, neqs, integrator, SD, integrator.p.SOL_VARS_TYPE)
        cfl_u = velomax*dt/Δs #Advective CFL
        cfl_c = c*dt/Δs       #Acoustic CFL

        Δs2      = Δs*Δs
        μ        = maximum(visc)
        λ        = 2.0 #free parameter
        cfl_visc = dt*λ*μ/Δs2 #Viscous CFL

        
        println_rank(" #  Advective CFL: ", cfl_u;    msg_rank = rank)
        println_rank(" #  Acoustic  CFL: ", cfl_c;    msg_rank = rank)
        println_rank(" #  Viscous   CFL: ", cfl_visc; msg_rank = rank) #, suppress = mesh.msg_suppress)
    end
end
