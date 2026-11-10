#----------------------------------------------------------------------------------------------
# TEMPORARY boundary/periodicity diagnostics, active only with inputs[:ldebug_bc] = true.
#
# To remove: delete this file, its include in src/Jexpresso.jl, the call in rhs.jl
# (search "debug_bc_rhs"), and the :ldebug_bc default in src/io/mod_inputs.jl.
#
# Called at the end of _build_rhs! (RHS is then the tendency, after global DSS and mass division).
# Prints, per variable:
#   * max |tendency| on: x-periodic face nodes, y-periodic face nodes, wall (non-periodic) face
#     nodes, and interior nodes, plus the location of the global max;
#   * max mismatch of the tendency between local copies of the same global node (periodic images
#     on this rank), which must be zero after a correct assembly;
#   * the domain integral of rho' (mass perturbation), counted once per global node.
#----------------------------------------------------------------------------------------------

const _DEBUG_BC_CALLS   = Ref(0)
_DBG_WORST = nothing
const _DEBUG_BC_CLASSES = Ref{Any}(nothing)

function _debug_bc_classify(mesh)
    npoin = mesh.npoin
    cls   = zeros(Int8, npoin)            # 0 interior, 1 periodicx, 2 periodicy, 3 periodicz, 4 wall
    for iface = 1:size(mesh.bdy_face_type, 1)
        tag = mesh.bdy_face_type[iface]
        c = tag in ("periodicx", "periodic1") ? Int8(1) :
            tag in ("periodicy", "periodic3") ? Int8(2) :
            tag in ("periodicz", "periodic2") ? Int8(3) : Int8(4)
        for i = 1:mesh.ngl, j = 1:mesh.ngl
            ip = mesh.poin_in_bdy_face[iface, i, j]
            # walls win over periodic for edge nodes, so pure periodic-face nodes stay separate
            cls[ip] = max(cls[ip], c)
        end
    end
    # local copies sharing a global id (periodic images on the same rank)
    groups = Dict{Int, Vector{Int}}()
    for ip = 1:npoin
        push!(get!(groups, Int(mesh.ip2gip[ip]), Int[]), ip)
    end
    dup = [v for v in values(groups) if length(v) > 1]
    return (cls = cls, dup = dup)
end

function debug_bc_rhs(params, time)
    (params.inputs[:ldebug_bc] || get(ENV, "JEXPRESSO_DEBUG_BC", "0") == "1") || return nothing
    _DEBUG_BC_CALLS[] += 1
    ncall = _DEBUG_BC_CALLS[]
    every = get(params.inputs, :ldebug_bc_every, 50)
    (ncall <= 5 || ncall % every == 0) || return nothing

    mesh = params.mesh
    if _DEBUG_BC_CLASSES[] === nothing
        _DEBUG_BC_CLASSES[] = _debug_bc_classify(mesh)
    end
    C    = _DEBUG_BC_CLASSES[]
    RHS  = params.RHS
    neqs = params.neqs
    comm = mesh.parts.comm
    rank = MPI.Comm_rank(comm)

    names = ("interior", "per-x", "per-y", "per-z", "wall")
    out = IOBuffer()
    println(out, " [ldebug_bc] call $ncall  t = $time  (dup groups on rank $rank: $(length(C.dup)))")
    vm = MPI.Allreduce([maximum(abs.(@view params.uaux[:, k])) for k in 2:4], MPI.MAX, comm)
    println(out, "   max|ρu'|,|ρv'|,|ρw'| = $(round.(vm; sigdigits = 3))")
    l2 = MPI.Allreduce([sum(abs2, @view params.uaux[:, k]) for k in 1:neqs], MPI.SUM, comm)
    println(out, "   sum(u^2) per var = $(round.(l2; sigdigits = 12))")
    # wall nodes: normal / tangential total momentum w.r.t. each wall face normal
    let mesh = params.mesh, mt = params.metrics, qe = params.qp.qe
        lpert = params.SOL_VARS_TYPE == PERT()
        mn = 0.0; mtg = 0.0
        for iface = 1:size(mesh.bdy_face_type, 1)
            tag = mesh.bdy_face_type[iface]
            tag in ("periodicx","periodic1","periodicy","periodic3","periodicz","periodic2") && continue
            for i = 1:mesh.ngl, j = 1:mesh.ngl
                ip = mesh.poin_in_bdy_face[iface, i, j]
                n1 = mt.nx[iface,i,j]; n2 = mt.ny[iface,i,j]; n3 = mt.nz[iface,i,j]
                m1 = params.uaux[ip,2] + (lpert ? qe[ip,2] : 0.0)
                m2 = params.uaux[ip,3] + (lpert ? qe[ip,3] : 0.0)
                m3 = params.uaux[ip,4] + (lpert ? qe[ip,4] : 0.0)
                d = m1*n1 + m2*n2 + m3*n3
                mn = max(mn, abs(d)); mtg = max(mtg, sqrt((m1-d*n1)^2 + (m2-d*n2)^2 + (m3-d*n3)^2))
            end
        end
        mn = MPI.Allreduce(mn, MPI.MAX, comm); mtg = MPI.Allreduce(mtg, MPI.MAX, comm)
        println(out, "   wall nodes: max|m.n| = $mn   max|m_t| = $mtg")
        if ncall == 1 && hasproperty(params, :bc_cons) && params.bc_cons !== nothing
            bc = params.bc_cons
            cnt = zeros(Int, 4)
            for kk = 1:length(bc.np_ip)
                P = @view bc.np_P[:, kk]
                r = round(Int, 3 - (P[1] + P[5] + P[9]))
                cnt[r + 1] += 1
            end
            println(out, "   no-penetration nodes by number of normals removed (0,1,2,3): $cnt")
            # deviation of node normals from the face-centre normal (planar faces): snapping check
            dev = 0.0; cm = (mesh.ngl + 1) ÷ 2
            for iface = 1:size(mesh.bdy_face_type, 1)
                tag = mesh.bdy_face_type[iface]
                tag in ("periodicx","periodic1","periodicy","periodic3","periodicz","periodic2") && continue
                c1 = mt.nx[iface,cm,cm]; c2 = mt.ny[iface,cm,cm]; c3 = mt.nz[iface,cm,cm]
                for i = 1:mesh.ngl, j = 1:mesh.ngl
                    d = abs(mt.nx[iface,i,j]-c1) + abs(mt.ny[iface,i,j]-c2) + abs(mt.nz[iface,i,j]-c3)
                    if d > dev
                        dev = d
                        global _DBG_WORST = (tag, i, j, (mt.nx[iface,i,j], mt.ny[iface,i,j], mt.nz[iface,i,j]), (c1,c2,c3),
                                             Tuple(mesh.coords[mesh.poin_in_bdy_face[iface,i,j], :]))
                    end
                end
            end
            println(out, "   max deviation of face-node normals from face-centre normal: $dev")
            dev > 1e-6 && println(out, "   worst: $(_DBG_WORST)")
        end
    end
    for ieq = 1:neqs
        m = zeros(5)
        for ip = 1:mesh.npoin
            c = C.cls[ip] + 1
            a = abs(RHS[ip, ieq])
            a > m[c] && (m[c] = a)
        end
        mism = 0.0
        for g in C.dup
            r0 = RHS[g[1], ieq]
            for ip in g
                mism = max(mism, abs(RHS[ip, ieq] - r0))
            end
        end
        m    = MPI.Allreduce(m, MPI.MAX, comm)
        mism = MPI.Allreduce(mism, MPI.MAX, comm)
        print(out, "   eq $ieq:")
        for c = 1:5
            print(out, "  ", names[c], " ", round(m[c]; sigdigits = 3))
        end
        println(out, "  | dup-mismatch ", round(mism; sigdigits = 3))
    end

    # mass perturbation integral, one contribution per global node (owned copy, first local image)
    Minv = params.Minv
    seen = Set{Int}()
    mass = 0.0
    for ip = 1:mesh.npoin
        g = Int(mesh.ip2gip[ip])
        (mesh.gip2owner[ip] == rank && !(g in seen)) || continue
        push!(seen, g)
        mass += params.uaux[ip, 1] / Minv[ip]
    end
    mass = MPI.Allreduce(mass, MPI.SUM, comm)
    println(out, "   ∫ρ' dV = $mass")

    # location of max |RHS u| and horizontal variation of p' and ρ' along z levels (first calls only)
    if ncall <= 1
        x = mesh.x; y = mesh.y; z = mesh.z
        _, imax = findmax(abs.(@view RHS[:, 2]))
        println(out, "   max|RHS u| at (x,y,z) = ($(x[imax]), $(y[imax]), $(z[imax])), class $(names[C.cls[imax]+1])")
        pp = params.uaux[:, end] .- params.qp.qe[:, end]
        rp = params.uaux[:, 1]
        lev = Dict{Float64, Vector{Int}}()
        for ip = 1:mesh.npoin
            push!(get!(lev, round(z[ip]; digits = 6), Int[]), ip)
        end
        zspread = maximum(maximum(z[v]) - minimum(z[v]) for v in values(lev))
        dp = maximum(maximum(pp[v]) - minimum(pp[v]) for v in values(lev))
        dr = maximum(maximum(rp[v]) - minimum(rp[v]) for v in values(lev))
        println(out, "   levels $(length(lev)), max z-spread in a level $zspread m")
        println(out, "   max|mesh.z - coords[:,3]| = $(maximum(abs.(z .- mesh.coords[:, 3])))   max|mesh.x - coords[:,1]| = $(maximum(abs.(x .- mesh.coords[:, 1])))")
        println(out, "   max|p'| $(maximum(abs.(pp)))  max horizontal spread of p' $dp")
        println(out, "   max|ρ'| $(maximum(abs.(rp)))  max horizontal spread of ρ' $dr")
        if isdefined(@__MODULE__, :model_eos_1pt)
            PC = PhysicalConst{Float64}(); MC = MicrophysicalConst{Float64}()
            qe = params.qp.qe
            d1 = 0.0; d2 = 0.0
            for ip = 1:mesh.npoin
                ρ = qe[ip, 1]
                Pm = Base.invokelatest(model_eos_1pt, ρ, qe[ip, 5] / ρ, qe[ip, 6] / ρ, z[ip], MC, PC)[1]
                d1 = max(d1, abs(Pm - qe[ip, end])); d2 = max(d2, abs(Pm - params.uaux[ip, end]))
            end
            println(out, "   max|P_1pt(qe) - qe[:,end]| = $d1   max|P_1pt(qe) - uaux[:,end]| = $d2")
        end
        println(out, "   max|P_ref| horizontal spread $(maximum(maximum(params.qp.qe[v, end]) - minimum(params.qp.qe[v, end]) for v in values(lev)))")
    end
    rank == 0 && print(String(take!(out)))
    return nothing
end
