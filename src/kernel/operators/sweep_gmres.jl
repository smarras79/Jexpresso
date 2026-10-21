"""
Estimates, per (iel,d), a scalar attenuation factor ρ ≈ exp(-κ·s_min),
where s_min is the SHORTEST straight-line distance for direction (iel,d)
to complete one period along WHICHEVER periodic axis it's most aligned
with (the dominant source of slow convergence). Used as a cheap,
analytic LEFT preconditioner: scale each direction's residual sub-vector
by 1/(1-ρ) before GMRES orthogonalizes against it.
"""
function compute_direction_preconditioner_scales(mesh, κ_mean::Float64, periodic_axes::Vector{Int})
    scales = Dict{Tuple{Int,Int}, Float64}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ω = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
            s_min = Inf
            for axis in periodic_axes
                Ωa = abs(Ω[axis])
                Ωa < 1e-6 && continue
                L_axis = axis==1 ? (mesh.xmax-mesh.xmin) : (axis==2 ? (mesh.ymax-mesh.ymin) : (mesh.zmax-mesh.zmin))
                s_min = min(s_min, L_axis/Ωa)
            end
            ρ = isfinite(s_min) ? exp(-κ_mean*s_min) : 0.0   # NOT periodic-aligned at all -> no slow mode
            scales[(iel,d)] = 1.0 / max(1.0 - ρ, 1e-3)        # capped to avoid dividing by ~0
        end
    end
    return scales
end

function apply_direction_precond!(flat::Vector{Float64}, scales_flat::Vector{Float64})
    flat .*= scales_flat
end


"""
Core restarted GMRES(restart) on an EXPLICIT matvec/RHS pair — no
knowledge of apply_G or preconditioning, so both plain and
preconditioned callers share this exact implementation.
"""
function gmres_core(matvec::Function, b::Vector{Float64}, x0::Vector{Float64}, comm;
                     restart::Int=20, max_restarts::Int=10, tol::Float64=1e-5, verbose::Bool=true)
    rank = MPI.Comm_rank(comm)
    n = length(x0)
    gdot(a,c) = MPI.Allreduce(dot(a,c), MPI.SUM, comm)
    gnorm(a) = sqrt(gdot(a,a))

    x = copy(x0)
    for outer = 1:max_restarts
        r0 = b .- matvec(x)
        β = gnorm(r0)
        rank==0 && verbose && println("  GMRES outer $outer: ||r0||=$β")
        β < tol && (rank==0 && println("GMRES converged (outer residual)"); return x)

        V = Vector{Vector{Float64}}(undef, restart+1)
        V[1] = r0 ./ β
        H = zeros(restart+1, restart)
        cs = zeros(restart); sn = zeros(restart)
        g = zeros(restart+1); g[1] = β

        m_used = restart
        for j = 1:restart
            w = matvec(V[j])
            for i = 1:j
                H[i,j] = gdot(w, V[i])
                w .-= H[i,j] .* V[i]
            end
            H[j+1,j] = gnorm(w)
            if H[j+1,j] < 1e-14
                m_used = j; break
            end
            V[j+1] = w ./ H[j+1,j]

            for i = 1:j-1
                temp      =  cs[i]*H[i,j] + sn[i]*H[i+1,j]
                H[i+1,j]  = -sn[i]*H[i,j] + cs[i]*H[i+1,j]
                H[i,j]    = temp
            end
            δ = sqrt(H[j,j]^2 + H[j+1,j]^2)
            cs[j] = H[j,j]/δ; sn[j] = H[j+1,j]/δ
            H[j,j] = δ; H[j+1,j] = 0.0
            temp = cs[j]*g[j]; g[j+1] = -sn[j]*g[j]; g[j] = temp

            rank==0 && verbose && println("    GMRES inner $j: residual=$(abs(g[j+1]))")
            abs(g[j+1]) < tol && (m_used = j; break)
        end

        y = zeros(m_used)
        for i = m_used:-1:1
            y[i] = (g[i] - sum(H[i,k]*y[k] for k = i+1:m_used; init=0.0)) / H[i,i]
        end
        for i = 1:m_used
            x .+= y[i] .* V[i]
        end
    end
    rank==0 && println("GMRES: reached max_restarts without meeting tol")
    return x
end

function gmres_fixed_point_solve(apply_G::Function, x0::Vector{Float64}, comm;
                                  restart::Int=20, max_restarts::Int=10, tol::Float64=1e-5, verbose::Bool=true)
    b = apply_G(zeros(length(x0)))
    matvec(v) = v .- (apply_G(v) .- b)
    return gmres_core(matvec, b, x0, comm; restart=restart, max_restarts=max_restarts, tol=tol, verbose=verbose)
end

function gmres_fixed_point_solve_precond(apply_G::Function, x0::Vector{Float64}, scales_flat::Vector{Float64}, comm;
                                          restart::Int=20, max_restarts::Int=10, tol::Float64=1e-5, verbose::Bool=true)
    b = apply_G(zeros(length(x0)))
    matvec(v) = scales_flat .* (v .- (apply_G(v) .- b))
    b_precond = scales_flat .* b
    return gmres_core(matvec, b_precond, x0, comm; restart=restart, max_restarts=max_restarts, tol=tol, verbose=verbose)
end

function flatten_scales(scales::Dict{Tuple{Int,Int},Float64}, order)
    len_each_lookup = Dict(k => length(v) for (k,v) in Dict{Tuple{Int,Int},Vector{Float64}}())  # placeholder, see below
    # Each (iel,d) contributes len_each IDENTICAL scale values (the scale
    # applies uniformly across the ngl^3 local nodes of that direction).
    error("needs len_each — see call site")
end

# Simpler, direct version used at the call site instead:
function build_scales_flat(scales::Dict{Tuple{Int,Int},Float64}, order, len_each::Int)
    vcat([fill(scales[k], len_each) for k in order]...)
end

"""
    detect_periodic_axes(mesh)

Scans mesh.bdy_face_type for "periodicx"/"periodicy"/"periodicz" tags and
returns the set of axes (1=x, 2=y, 3=z) that are actually periodic on
this mesh. Reused wherever an axis-aware periodic computation is needed
(currently just the GMRES direction preconditioner).
"""
function detect_periodic_axes(mesh)
    axes = Set{Int}()
    axis_char_to_int = Dict('x'=>1, 'y'=>2, 'z'=>3)
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        tag = lowercase(string(mesh.bdy_face_type[iface_bdy]))
        occursin("periodic", tag) || continue
        for (ch, ax) in axis_char_to_int
            occursin("periodic$ch", tag) && push!(axes, ax)
        end
    end
    return collect(axes)
end

function jacobi_with_extrapolation(apply_G::Function, x0::Vector{Float64}, comm;
                                    len_each::Int=0, max_iter::Int=300, tol::Float64=1e-8,
                                    mode::Symbol=:global,
                                    cos_thresh::Float64=0.9, max_factor::Float64=50.0,
                                    jump_scale::Float64=0.5, accept_ratio::Float64=1.0,
                                    reject_backoff::Int=2, max_backoff::Int=16,
                                    verbose::Bool=true)
    rank = MPI.Comm_rank(comm)
    gdot(a,b) = MPI.Allreduce(dot(a,b), MPI.SUM, comm)
    gmax(a) = MPI.Allreduce(isempty(a) ? 0.0 : maximum(abs.(a)), MPI.MAX, comm)

    function build_candidate(x, d, dprev)
        if mode == :global
            dd = gdot(d, dprev); pp = gdot(dprev, dprev); nn = gdot(d, d)
            λ = dd / max(pp, 1e-300)
            cs = dd / sqrt(max(nn*pp, 1e-300))
            rank == 0 && verbose && println("    global fit: cos=$cs, λ=$λ")
            (cs > cos_thresh && 0 < λ < 1) || return nothing
            return x .+ min(jump_scale*λ/(1-λ), max_factor) .* d
        else  # :block — one λ per (iel,d) block of len_each entries (all local, no communication)
            nb = length(d) ÷ len_each
            cand = copy(x); n_active = 0; fsum = 0.0
            for b = 1:nb
                r = (b-1)*len_each+1 : b*len_each
                dv = view(d, r); pv = view(dprev, r)
                pp = dot(pv, pv); pp < 1e-300 && continue
                dd = dot(dv, pv); nn = dot(dv, dv)
                λ = dd/pp; cs = dd/sqrt(max(nn*pp, 1e-300))
                (cs > cos_thresh && 0 < λ < 1) || continue
                f = min(jump_scale*λ/(1-λ), max_factor)
                cand[r] .+= f .* dv
                n_active += 1; fsum += f
            end
            tot_active = MPI.Allreduce(n_active, MPI.SUM, comm)
            rank == 0 && verbose && println("    block fit: $tot_active blocks active (rank 0 mean factor $(n_active>0 ? fsum/n_active : 0.0))")
            return tot_active > 0 ? cand : nothing   # tot_active is identical on every rank
        end
    end

    x = copy(x0)
    dprev = nothing; dlast = nothing
    n_sweeps = 0; n_acc = 0; n_rej = 0
    backoff = 0; cur_backoff = reject_backoff

    while n_sweeps < max_iter
        g = apply_G(x); n_sweeps += 1
        d = g .- x
        dmax = gmax(d)
        x = g
        dprev = dlast; dlast = d
        rank == 0 && verbose && println("  JX sweep $n_sweeps: max_delta=$dmax")
        if dmax < tol
            rank == 0 && println("JX converged: $n_sweeps sweeps, $n_acc accepted, $n_rej rejected")
            return x
        end
        backoff > 0 && (backoff -= 1)

        while dprev !== nothing && backoff == 0 && n_sweeps < max_iter
            cand = build_candidate(x, dlast, dprev)
            cand === nothing && break
            gc = apply_G(cand); n_sweeps += 1
            dc = gc .- cand
            dcmax = gmax(dc); dlastmax = gmax(dlast)
            rank == 0 && verbose && println("    jump attempt: delta at candidate / previous delta = $(dcmax/max(dlastmax,1e-300))")
            if dcmax < accept_ratio * dlastmax
                n_acc += 1; cur_backoff = reject_backoff
                x = gc
                if dcmax < tol
                    rank == 0 && println("JX converged: $n_sweeps sweeps, $n_acc accepted, $n_rej rejected")
                    return x
                end
                dlast = dc; dprev = nothing
                break
            else
                n_rej += 1
                backoff = cur_backoff
                cur_backoff = min(2*cur_backoff, max_backoff)
                break
            end
        end
    end
    rank == 0 && println("JX: sweep budget reached ($n_sweeps sweeps, $n_acc accepted, $n_rej rejected)")
    return x
end
