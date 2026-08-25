function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        println(" Initialize fields for 3D CompEuler with θ equation ........................ ")
    end

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρw", "ρθ"]
    qoutvars = ["ρ", "u", "v", "w", "θ", "p"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # A SPHERE, not a cylinder. CompEuler/rtb_imex measures r in the x-z plane
    # only, which makes its bubble a cylinder along y and its solution
    # y-invariant -- deliberately, because that case exists to reproduce the 2D
    # CompEuler/theta as a slab.
    #
    # This one is genuinely three-dimensional: r is the full 3D radius, so the
    # bubble is a ball and the flow it drives has structure in all three
    # directions. That matters for what this case is FOR. The Schur reduction
    # eliminates rho_u and rho_v along with everything else, and a run where
    # both are identically zero would exercise the elimination on a degenerate
    # state and report a saving that a real 3D flow might not repeat.
    #
    # Everything else is rtb_imex's bubble unchanged -- 2 K amplitude, 2000 m
    # radius, linear taper theta_c(1 - r/r0), 300 K background, centred at
    # z = 2500 m -- so the physics is one already exercised in this repository
    # and only the dimensionality is new.
    #
    # NOTE ON WHAT THE AMPLITUDE DOES NOT CONTROL: how diffused the result
    # looks is set by sqrt(mu*t)/r0, not by theta_c. At mu = 125 and t = 1000 s
    # on r0 = 2000 m that is 0.18, matching rtb_imex.
    #-------------------------------------------------------------------------
    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        if inputs[:lrestart] == true
            #
            # READ RESTART HDF5:
            #
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))

            for ip=1:mesh.npoin
                ρ  = q.qn[ip,1]
                ρθ = q.qn[ip,5]
                θ  = ρθ/ρ
                P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
                q.qn[ip,end] = P

                ρe  = q.qe[ip,1]
                ρθe = q.qe[ip,5]
                θe  = ρθe/ρ
                Pe = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
                q.qe[ip,end] = Pe
            end

        else
            #
            # INITIAL STATE from scratch:
            #
            # The centre is taken from the MESH rather than hard-coded, and
            # through Allreduce rather than from the local extrema: on more
            # than one rank `maximum(mesh.x)` is the local maximum, so a
            # hard-coded or rank-local centre would put a different bubble on
            # every rank and the answer would depend on the partition.
            max_x = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
            min_x = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
            max_y = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm)
            min_y = MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)

            xc = (max_x + min_x)/2
            yc = (max_y + min_y)/2
            zc = 2500.0 #m
            r0 = 2000.0 #m

            θref = 300.0 #K
            θc   =   2.0 #K
            for ip = 1:mesh.npoin

                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]

                r = sqrt( (x - xc)^2 + (y - yc)^2 + (z - zc)^2 )

                Δθ = 0.0 #K
                if r < r0
                    Δθ = θc*(1.0 - r/r0)
                end
                θ = θref + Δθ
                p    = PhysConst.pref*(1.0 - PhysConst.g*z/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
                pref = PhysConst.pref*(1.0 - PhysConst.g*z/(PhysConst.cp*θref))^(PhysConst.cpoverR)
                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
                ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

                u = 0.0
                v = 0.0
                w = 0.0

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*w - ρref*w
                    q.qn[ip,5] = ρ*θ - ρref*θref
                    q.qn[ip,end] = p

                    #Store initial background state for plotting and analysis of perturbations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = u
                    q.qe[ip,3] = v
                    q.qe[ip,4] = w
                    q.qe[ip,5] = ρref*θref
                    q.qe[ip,end] = pref
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*w
                    q.qn[ip,5] = ρ*θ
                    q.qn[ip,end] = p

                    #Store initial background state for plotting and analysis of perturbations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u
                    q.qe[ip,3] = ρref*v
                    q.qe[ip,4] = ρref*w
                    q.qe[ip,5] = ρref*θref
                    q.qe[ip,end] = pref
                end
            end
        end

        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,5] .= q.qn[:,5]./(q.qn[:,1] + q.qe[:,1])

                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
                q.qn[:,5] .= q.qn[:,5]./q.qn[:,1]

                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            end
        end

    else
        lpert = (inputs[:SOL_VARS_TYPE] == PERT())
        PhysConst = PhysicalConst{TFloat}()
        # Same Allreduce as the CPU branch, for the same reason.
        xc = TFloat((MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm) +
                     MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm))/2)
        yc = TFloat((MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm) +
                     MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm))/2)
        zc = TFloat(2500.0) #m
        rθ = TFloat(2000.0) #m

        θref = TFloat(300.0) #K
        θc   =   TFloat(2.0) #K
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, mesh.z, xc, yc, rθ, zc, θref, θc, PhysConst, lpert; ndrange = (mesh.npoin))
    end
    if rank == 0
        println(" Initialize fields for 3D CompEuler with θ equation ........................ DONE ")
    end

    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, z, xc, yc, rθ, zc, θref, θc, PhysConst, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x = x[ip]
    y = y[ip]
    z = z[ip]
    # 3D radius, matching the CPU branch. A 2D radius here and a 3D one there
    # would give a different initial state on GPU and CPU, and nothing in the
    # output would say which one was run.
    r = sqrt( (x - xc)^2 + (y - yc)^2 + (z - zc)^2 )
    Δθ = T(0.0) #K

    if r < rθ
        Δθ = T(θc*(T(1.0) - r/rθ))
    end

    θ = θref + Δθ
    p    = PhysConst.pref*(T(1.0) - PhysConst.g*z/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
    pref = PhysConst.pref*(T(1.0) - PhysConst.g*z/(PhysConst.cp*θref))^(PhysConst.cpoverR)
    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

    u = T(0.0)
    v = T(0.0)
    w = T(0.0)

    if (lpert)
        qn[ip,1] = ρ - ρref
        qn[ip,2] = ρ*u - ρref*u
        qn[ip,3] = ρ*v - ρref*v
        qn[ip,4] = ρ*w - ρref*w
        qn[ip,5] = ρ*θ - ρref*θref
        qn[ip,end] = p
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*w
        qn[ip,5] = ρ*θ
        qn[ip,end] = p
    end

    qe[ip,1] = ρref
    qe[ip,2] = ρref*u
    qe[ip,3] = ρref*v
    qe[ip,4] = ρref*w
    qe[ip,5] = ρref*θref
    qe[ip,end] = pref
end

function user_get_adapt_flags!(adapt_flags, inputs, old_ad_lvl, q, qe, connijk, nelem, ngl)
    # AMR is off for this case -- IMEX3D refuses :ladapt, because adaptation
    # invalidates the column topology the preconditioner is built on. Kept so
    # the case has the same shape as its neighbours.
    ips         = KernelAbstractions.zeros(CPU(), TInt, ngl * ngl * ngl)
    tol         = 301.2
    max_level   = inputs[:amr_max_level]

    for iel = 1:nelem
        m = 1
        for i = 1:ngl, j = 1:ngl, k = 1:ngl
            ips[m] = connijk[iel, i, j, k]
            m += 1
        end
        theta = q[ips, 5] ./ q[ips, 1]
        if any(theta .> tol) && (old_ad_lvl[iel] < max_level)
            adapt_flags[iel] = refine_flag
        end
        if all(theta .< tol)
            adapt_flags[iel] = coarsen_flag
        end
    end
end
