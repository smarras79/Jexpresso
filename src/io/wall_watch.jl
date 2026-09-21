#
# wall_watch_report -- one line per call on the state of the MOST wall nodes.
#
# Collective (Allreduce), so every rank must call it. Enabled with
# JEXPRESSO_WALL_WATCH=<steps> in the environment; see TimeIntegrators.jl.
#
# Prints, for the wall node with the largest |u_h| anywhere in the domain:
# where it is, its own velocity, the velocity at the first node above it
# (ifirst_wall_node_index), rho, theta on both, mu_turb on both, and the
# u_star / tau that MOST would build from the node above. Then a summary over
# every wall node: how many are above `thresh` m/s, the mean |u_h| offset
# between the wall node and the node above, and the extrema on the layer
# above. Meant to be pasted, not parsed.
#
function wall_watch_report(params, u, t, step; io = stdout, thresh = 12.0)
    mesh   = params.mesh
    inputs = params.inputs
    comm   = MPI.COMM_WORLD
    rank   = MPI.Comm_rank(comm)
    neqs   = Int(params.neqs)

    u2uaux!(@view(params.uaux[:,:]), u, neqs, mesh.npoin)
    uaux   = params.uaux
    qe     = params.qp.qe
    lpert  = inputs[:SOL_VARS_TYPE] == PERT()
    ifw    = inputs[:ifirst_wall_node_index]::Int
    ngl    = Int(mesh.ngl)
    coords = mesh.coords
    have_μ = params.sgs !== nothing && hasproperty(params.sgs, :μ_turb) &&
             length(params.sgs.μ_turb) == mesh.npoin
    μt     = have_μ ? params.sgs.μ_turb : nothing
    κ      = PhysicalConst{Float64}().karman
    z0m    = 0.1        # what BCs.jl passes to CM_MOST! in the dry branch

    @inline function state(ip)
        if lpert
            ρ = uaux[ip,1] + qe[ip,1]
            return ρ, (uaux[ip,2]+qe[ip,2])/ρ, (uaux[ip,3]+qe[ip,3])/ρ,
                      (uaux[ip,4]+qe[ip,4])/ρ, (uaux[ip,5]+qe[ip,5])/ρ
        else
            ρ = uaux[ip,1]
            return ρ, uaux[ip,2]/ρ, uaux[ip,3]/ρ, uaux[ip,4]/ρ, uaux[ip,5]/ρ
        end
    end

    best   = -Inf
    rec    = zeros(Float64, 14)     # x y u v u1 v1 w1 rho th th1 mu mu1 z1 rho1
    n_wall = 0; n_run = 0; sum_off = 0.0
    max_u1 = 0.0; max_w1 = 0.0; min_ρθ = Inf

    for iface = 1:mesh.nfaces_bdy
        mesh.bdy_face_type[iface] == "MOST" || continue
        e = mesh.bdy_face_in_elem[iface]
        for i = 1:ngl, j = 1:ngl
            ip  = mesh.poin_in_bdy_face[iface,i,j]
            ip1 = mesh.connijk[e,i,j,ifw]
            ρ,  uu,  vv,  ww,  th  = state(ip)
            ρ1, uu1, vv1, ww1, th1 = state(ip1)
            uh  = hypot(uu,  vv)
            uh1 = hypot(uu1, vv1)
            n_wall  += 1
            sum_off += abs(uh - uh1)
            uh > thresh && (n_run += 1)
            max_u1 = max(max_u1, uh1)
            max_w1 = max(max_w1, abs(ww1))
            min_ρθ = min(min_ρθ, ρ*th, ρ1*th1)
            if uh > best
                best = uh
                rec[1]  = coords[1,ip];  rec[2]  = coords[2,ip]
                rec[3]  = uu;            rec[4]  = vv
                rec[5]  = uu1;           rec[6]  = vv1;   rec[7] = ww1
                rec[8]  = ρ;             rec[9]  = th;    rec[10] = th1
                rec[11] = have_μ ? μt[ip] : NaN
                rec[12] = have_μ ? μt[ip1] : NaN
                rec[13] = coords[3,ip1] - coords[3,ip]
                rec[14] = ρ1
            end
        end
    end

    gbest  = MPI.Allreduce(best, MPI.MAX, comm)
    owner  = MPI.Allreduce(best == gbest ? rank : typemax(Int), MPI.MIN, comm)
    g_n    = MPI.Allreduce(n_wall, +, comm)
    g_run  = MPI.Allreduce(n_run,  +, comm)
    g_off  = MPI.Allreduce(sum_off, +, comm)
    g_u1   = MPI.Allreduce(max_u1, MPI.MAX, comm)
    g_w1   = MPI.Allreduce(max_w1, MPI.MAX, comm)
    g_ρθ   = MPI.Allreduce(min_ρθ, MPI.MIN, comm)
    # Peak resident set size across ranks, in GB. mpirun-launched ranks are
    # invisible to sacct, and job 1300883 died of std::bad_alloc at t = 9540
    # with no record of how memory had grown; this is the record.
    g_rss  = MPI.Allreduce(Float64(Sys.maxrss()) / 2^30, MPI.MAX, comm)

    if rank == owner && isfinite(gbest)
        uh1   = hypot(rec[5], rec[6])
        ustar = κ * uh1 / log(rec[13] / z0m)
        τ     = rec[14] * ustar^2
        @printf(io, " # wall-watch t=%.1f step=%d rank=%d | WALL max|uh|=%.2f at (x,y)=(%.0f,%.0f) u,v=(%.2f,%.2f) rho=%.3f th=%.2f mu_t=%.3g",
                Float64(t), step, rank, gbest, rec[1], rec[2], rec[3], rec[4], rec[8], rec[9], rec[11])
        @printf(io, " | NODE2 (z1=%.2f) u,v,w=(%.2f,%.2f,%.2f) th=%.2f mu_t=%.3g | MOST u*=%.3f tau=%.3f\n",
                rec[13], rec[5], rec[6], rec[7], rec[10], rec[12], ustar, τ)
        flush(io)
    end
    MPI.Barrier(comm)
    if rank == 0
        @printf(io, " # wall-watch t=%.1f summary | wall nodes=%d, |uh|>%.0f: %d | mean|uh_wall-uh_node2|=%.3f | node2 layer max|uh|=%.2f max|w|=%.2f | min(rho*theta)=%.1f | max RSS=%.2f GB\n",
                Float64(t), g_n, thresh, g_run, g_off / max(g_n,1), g_u1, g_w1, g_ρθ, g_rss)
        flush(io)
    end
    return nothing
end
