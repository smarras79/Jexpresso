using Printf, Dates

const PROFILE_CODE_TAG = "prof-v18"                                        # bump when the profiling code changes
const PROFILE_SRC_HASH = try string(hash(read(@__FILE__)); base=16) catch; "unavailable" end

Base.@kwdef struct RunConfig
    label::String = "run"
    gc_mode::Symbol = :auto
    gc_every::Int = 0
    local_passes::Bool = true
    measure_skew::Bool = false
    coupling::Symbol = :isotropic
    n_outer::Int = 2
    warmup_outer::Int = 1
    frac::Float64 = 0.01
    seed::Int = 42
    max_e_ext::Int = 4
    κ_val::Float64 = 0.001
    σ_val::Float64 = 0.0005
    face_plans::Bool = true
    gs_mode::Symbol = :cached
    spatial::Bool = true
    angular::Bool = true
    lu_cache_mb::Float64 = 0.0
    local_plans::Bool = true
    flat_faces::Bool = true
    angular_plans::Bool = true
    upwind_order::Bool = true
    incremental_pack::Bool = true
    round_cap::Int = 0
    entry_barrier::Bool = false
    gs_inplace::Bool = true
    gs_gc_off::Bool = true
    cache_min_pairs::Int = 30000
end

mutable struct SweepLog
    passes::Vector{Int}; solved::Vector{Int}
    t_compute::Vector{Float64}; t_exchange::Vector{Float64}; t_skew::Vector{Float64}; t_reduce::Vector{Float64}
    sent::Vector{Int}; valid::Vector{Int}; fresh::Vector{Int}; t_empty::Vector{Float64}
    gc_compute::Vector{Float64}; gc_exchange::Vector{Float64}          # NEW
end
SweepLog() = SweepLog(Int[], Int[], Float64[], Float64[], Float64[], Float64[], Int[], Int[], Int[], Float64[],
                      Float64[], Float64[])


mutable struct SweepProfile
    t::Dict{Symbol,Float64}
    n::Dict{Symbol,Int}
    sweeps::Vector{SweepLog}
end
SweepProfile() = SweepProfile(Dict{Symbol,Float64}(), Dict{Symbol,Int}(), SweepLog[])

@inline ptic() = time_ns()
@inline function padd!(p::SweepProfile, key::Symbol, dt::Float64)
    p.t[key] = get(p.t, key, 0.0) + dt; p.n[key] = get(p.n, key, 0) + 1; return nothing
end
@inline ptoc!(p::SweepProfile, key::Symbol, t0::UInt64) = padd!(p, key, (time_ns() - t0) * 1e-9)

const FEATURE_NAMES = ["nelem","pairs","exch_entries","xrank_fine","xrank_coarse","xrank_conf",
                       "xrank_periodic","ncf_child","ncf_parent","n_neighbor_ranks"]

# Fixed key list, so every rank enters the same reductions in the same order.
const PROF_SECTIONS = [:sweep_total, :gs_total,
    :attempt_solved, :attempt_wasted, :empty_pass,                       # empty_pass is a subset of attempt_wasted
    :classify, :face_upwind, :trial_eval, :alloc_A, :vol_matrix, :source_term, :flux_terms, :solve,
    :ex_pack, :ex_mpi, :ex_unpack, :skew_exchange, :skew_reduce, :allreduce, :bookkeeping, :periodic_exch,
    :outer_src, :outer_update, :sync_gc]
const BRANCH_LABELS = [:boundary, :periodic_same, :periodic_cross, :same_rank, :ncf_child, :ncf_parent,
                       :xrank_fine, :xrank_coarse, :xrank_angular, :xrank_ordinary]
const BRANCH_READY   = Dict(b => Symbol("fu_$(b)_ready")   for b in BRANCH_LABELS)
const BRANCH_BLOCKED = Dict(b => Symbol("fu_$(b)_blocked") for b in BRANCH_LABELS)
const GS_BRANCH_UP = Dict(b => Symbol("gsu_$(b)") for b in BRANCH_LABELS)      # top level, after BRANCH_LABELS
const GS_BRANCH_EV = Dict(b => Symbol("gse_$(b)") for b in BRANCH_LABELS)
for b in BRANCH_LABELS; push!(PROF_SECTIONS, GS_BRANCH_UP[b]); push!(PROF_SECTIONS, GS_BRANCH_EV[b]); end

for b in BRANCH_LABELS; push!(PROF_SECTIONS, BRANCH_READY[b]); push!(PROF_SECTIONS, BRANCH_BLOCKED[b]); end

push!(PROF_SECTIONS, :entry_wait)
append!(PROF_SECTIONS, [:run_fields, :run_init, :run_cache])
@inline gtic(p) = p === nothing ? UInt64(0) : time_ns()
@inline gtoc!(p, k::Symbol, t0::UInt64) = (p === nothing || ptoc!(p, k, t0); nothing)
append!(PROF_SECTIONS, [:gs_build, :gs_upwind, :gs_solve, :gs_exchange, :gs_reduce])
push!(PROF_SECTIONS, :entry_wait2)

function rank_features(setup)
    cnt(k) = count(v -> v == k, values(setup.neighbor_kind))
    return Float64[setup.mesh.nelem, sum(setup.n_dirs), setup.gid_exch.n_flat, cnt(:i_am_fine), cnt(:i_am_coarse),
                   cnt(:conforming), count(values(setup.neighbor_is_periodic)), length(setup.ncf_as_child),
                   length(setup.ncf_as_parent), length(unique(values(setup.neighbor_rank)))]
end
function pearson(a, b)
    da = a .- sum(a)/length(a); db = b .- sum(b)/length(b); d = sqrt(sum(da.^2)*sum(db.^2))
    return d > 0 ? sum(da .* db)/d : NaN
end


function per_gel_fingerprint_csv(setup, I_final, comm, path::String)
    mesh = setup.mesh; ngl = mesh.ngl; rank = MPI.Comm_rank(comm)
    cnt = zeros(Int, mesh.nelem, 5)      # cross-rank conforming, cross-rank NCF, cross-rank periodic, same-rank NCF, same-rank periodic
    for ((iel, fid), kind) in setup.neighbor_kind
        col = get(setup.neighbor_is_periodic, (iel, fid), false) ? 3 : (kind == :conforming ? 1 : 2)
        cnt[Int(iel), col] += 1
    end
    for (iel, _) in keys(setup.ncf_as_child);  cnt[Int(iel), 4] += 1; end
    for (iel, _) in keys(setup.ncf_as_parent); cnt[Int(iel), 4] += 1; end
    for ((iel, _), v) in setup.same_rank_is_periodic; v && (cnt[Int(iel), 5] += 1); end
    rows = Float64[]
    for iel = 1:mesh.nelem
        gel = Float64(mesh.el2gel[iel]); em = mesh.extra_mesh[iel]
        s = 0.0; c = 0.0; mx = 0.0
        for d = 1:em.extra_npoin
            v = I_final[(iel, d)]
            for p in eachindex(v)
                y = sin(0.7*gel + 1.3*d + 0.017*p) * v[p] - c
                t = s + y; c = (t - s) - y; s = t
                mx = max(mx, abs(v[p]))
            end
        end
        cx = 0.0; cy = 0.0; cz = 0.0; xlo = Inf; xhi = -Inf
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j, k]
            cx += mesh.x[ip]; cy += mesh.y[ip]; cz += mesh.z[ip]
            xlo = min(xlo, mesh.x[ip]); xhi = max(xhi, mesh.x[ip])
        end
        append!(rows, Float64[gel, rank, cx/ngl^3, cy/ngl^3, cz/ngl^3, xhi - xlo, s, mx,
                              cnt[iel,1], cnt[iel,2], cnt[iel,3], cnt[iel,4], cnt[iel,5]])
    end
    all = gather_all(rows, comm)                         # collective
    rank == 0 || return String[]
    M = permutedims(reshape(all, 13, :)); M = M[sortperm(M[:,1]), :]
    open(path, "w") do io
        println(io, "gel,rank,cx,cy,cz,xextent,wsum,maxabs,xrank_conf,xrank_ncf,xrank_periodic,same_ncf,same_periodic")
        for r in 1:size(M,1); println(io, join(M[r,:], ",")); end
    end
    return String[@sprintf("per-gel fingerprint: %d elements, sum of per-gel weighted sums %.15g (should equal weighted_sum above to ~1e-10)", size(M,1), sum(M[:,7]))]
end

function straggler_diagnosis(prof::SweepProfile, setup, comm, mats, ni, base::String)
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    tot(f) = sum(sum(getfield(L, f)) for L in prof.sweeps)
    bk = vcat([Symbol("fu_$(b)_$(s)") for b in BRANCH_LABELS for s in ("ready","blocked")],
          [BRANCH_EVAL[b] for b in BRANCH_EVAL_LABELS])
    v = vcat(rank_features(setup),
             Float64[tot(:t_compute), tot(:t_exchange), tot(:t_reduce), get(prof.t,:attempt_solved,0.0),
                     get(prof.t,:attempt_wasted,0.0), get(prof.t,:face_upwind,0.0)],
             Float64[get(prof.t,k,0.0) for k in bk], Float64[get(prof.n,k,0) for k in bk])
    K = length(v); nf = length(FEATURE_NAMES); nb = length(bk)
    M = reshape(gather_all(v, comm), K, np)                       # column k = rank k-1
    rank == 0 || return String[]

    med(x) = sort(x)[(end+1)÷2]
    comp = M[nf+1,:]; order = sortperm(comp; rev=true)[1:min(5,np)]; w = order[1]
    out = String["straggler diagnosis (per-rank totals over all measured sweeps):",
        "  correlation of per-rank compute time with: " *
            join(["$(FEATURE_NAMES[i])=$(round(pearson(comp, M[i,:]); digits=2))" for i in 1:nf], " ")]
    push!(out, "  5 busiest ranks by compute [rank node | " * join(FEATURE_NAMES[1:nf], ",") * " | compute attempt_wasted face_upwind s]:")
    for k in order
        push!(out, "    $(k-1) $(ni.node_of_rank[k]) | " * join(Int.(M[1:nf,k]), ",") *
                   @sprintf(" | %.3f %.3f %.3f", M[nf+1,k], M[nf+5,k], M[nf+6,k]))
    end
    push!(out, "  mean over ranks of the same three: " * @sprintf("%.3f %.3f %.3f", sum(M[nf+1,:])/np, sum(M[nf+5,:])/np, sum(M[nf+6,:])/np))
    push!(out, "  per-branch timers (fu_* = upwind construction, fe_* = closure evaluation): mean-rank s | busiest rank (s, calls, us/call):")
    for (j,k) in enumerate(bk)
        tm = M[nf+6+j,:]; cl = M[nf+6+nb+j,:]; sum(cl) == 0 && continue
        push!(out, @sprintf("    %-26s %9.3f | %9.3f %9d %9.1f", String(k), sum(tm)/np, tm[w], Int(cl[w]), 1e6*tm[w]/max(cl[w],1)))
    end
    m = mats[1]
    push!(out, "  busiest rank $(w-1), sweep 1, per round: round compute exchange reduce solved passes | median-rank compute")
    for r in 1:m.nr
        push!(out, @sprintf("    %3d %8.3f %8.3f %8.3f %6d %3d | %8.4f", r, m.C[r,w], m.X[r,w], m.R[r,w], Int(m.S[r,w]), Int(m.P[r,w]), med(m.C[r,:])))
    end
    names = vcat(FEATURE_NAMES, ["compute","exchange","reduce","attempt_solved","attempt_wasted","face_upwind"],
                 ["t_" * String(k) for k in bk], ["n_" * String(k) for k in bk])
    open(base * "_ranktable.csv", "w") do io
        println(io, "rank,node," * join(names, ","))
        for k in 1:np; println(io, "$(k-1),$(ni.node_of_rank[k])," * join(M[:,k], ",")); end
    end
    open(base * "_rankrounds.csv", "w") do io
        println(io, "sweep,round,rank,compute,exchange,reduce,solved,passes,gc_compute,gc_exchange")
        for (s, mm) in enumerate(mats), r in 1:mm.nr, k in 1:np
            println(io, join(Any[s, r, k-1, mm.C[r,k], mm.X[r,k], mm.R[r,k], Int(mm.S[r,k]), Int(mm.P[r,k]), mm.GC1[r,k], mm.GC2[r,k]], ","))
        end
    end
    return out
end


function face_branch(neighbor_elem, neighbor_rank, neighbor_kind, ncf_as_child, ncf_as_parent,
                     same_rank_is_periodic, neighbor_is_periodic, angular_deps_cr, iel, fid, d)
    get(same_rank_is_periodic, (iel,fid), false) && return :periodic_same
    get(neighbor_is_periodic,  (iel,fid), false) && return :periodic_cross
    neighbor_elem[iel][fid] != zero(TInt) && return :same_rank
    haskey(ncf_as_child,  (iel,fid)) && return :ncf_child
    haskey(ncf_as_parent, (iel,fid)) && return :ncf_parent
    if haskey(neighbor_rank, (iel,fid))
        k = neighbor_kind[(iel,fid)]
        k == :i_am_fine   && return :xrank_fine
        k == :i_am_coarse && return :xrank_coarse
        return haskey(angular_deps_cr, (iel,fid,d)) ? :xrank_angular : :xrank_ordinary
    end
    return :boundary
end



function config_from_env()
    e(k, d) = get(ENV, "PROFILE_" * k, d)
    return RunConfig(label=e("LABEL","run"), gc_mode=Symbol(e("GC_MODE","auto")), gc_every=parse(Int, e("GC_EVERY","0")),
        local_passes=parse(Bool, e("LOCAL_PASSES","true")), measure_skew=parse(Bool, e("MEASURE_SKEW","false")),
        coupling=Symbol(e("COUPLING","isotropic")), n_outer=parse(Int, e("N_OUTER","2")),
        warmup_outer=parse(Int, e("WARMUP_OUTER","1")), frac=parse(Float64, e("FRAC","0.01")),
        seed=parse(Int, e("SEED","42")), κ_val=parse(Float64, e("KAPPA","0.001")), σ_val=parse(Float64, e("SIGMA","0.0005")),
	face_plans=parse(Bool, e("FACE_PLANS","true")), gs_mode=Symbol(e("GS_MODE","cached")),
	upwind_order=parse(Bool, e("UPWIND_ORDER","true")), angular_plans=parse(Bool, e("ANGULAR_PLANS","true")),
	local_plans=parse(Bool, e("LOCAL_PLANS","true")), flat_faces=parse(Bool, e("FLAT_FACES","true")),
	cache_min_pairs=parse(Bool, e("CACHE_MIN_PAIRS","true")), gs_inplace=parse(Bool, e("GS_INPLACE","true")), gs_gc_off=parse(Bool, e("GS_GC_OFF","true")),
	incremental_pack=parse(Bool, e("INCREMENTAL_PACK","true")), round_cap=parse(Int, e("ROUND_CAP","0")), entry_barrier=parse(Bool, e("ENTRY_BARRIER","false")),
	spatial=parse(Bool, e("SPATIAL","true")), angular=parse(Bool, e("ANGULAR","true")), lu_cache_mb=parse(Float64, e("LU_CACHE_MB","0")))
end

mutable struct RunReport
    io::Union{IO,Nothing}
end
function emit(r::RunReport, s::AbstractString)
    r.io === nothing && return
    println(r.io, s); println(stdout, s); flush(r.io); flush(stdout)
end

const PROFILE_ENV_VARS = ["SLURM_NODELIST", "JULIA_NUM_THREADS","JULIA_NUM_GC_THREADS","JULIA_EXCLUSIVE","OMP_NUM_THREADS","MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS","I_MPI_PIN","I_MPI_PIN_DOMAIN","I_MPI_FABRICS","SLURM_CPU_BIND","SLURM_NNODES","SLURM_NTASKS","SLURM_JOB_ID"]

function gc_snapshot()
    g = Base.gc_num()
    return Float64[Base.gc_time_ns()*1e-9, hasproperty(g, :full_sweep) ? Float64(g.full_sweep) : NaN, 0.0, Float64(g.total_allocd)]
end

function gather_sweep_matrices(prof::SweepProfile, s::Int, comm)          # collective; column k = rank k-1
    np = MPI.Comm_size(comm); L = prof.sweeps[s]; nr = length(L.t_compute)
    has_gc = hasfield(SweepLog, :gc_compute) && length(L.gc_compute) == nr && length(L.gc_exchange) == nr
    mk(v) = reshape(gather_all(Float64.(v), comm), nr, np)
    return (nr=nr, has_gc=has_gc, C=mk(L.t_compute), X=mk(L.t_exchange), S=mk(L.solved), P=mk(L.passes), R=mk(L.t_reduce),
            GC1=mk(has_gc ? L.gc_compute : zeros(nr)), GC2=mk(has_gc ? L.gc_exchange : zeros(nr)))
end

function straggler_lines(m, node_of_rank; top::Int=5)
    nr = m.nr; np = size(m.C, 2)
    med(v) = sort(v)[(end+1)÷2]
    out = String[]
    m.has_gc || push!(out, "  NOTE: per-round GC arrays missing; all gc values below are zeros")
    function worst(T, G, label, with_solved)
        order = sortperm(vec(maximum(T; dims=2)); rev=true)[1:min(top, nr)]
        push!(out, "  $label, worst $(length(order)) rounds (slowest rank; gc = its GC time in the same phase):")
        n_gc = 0
        for r in order
            k = argmax(T[r,:]); dom = G[r,k] > 0.5*T[r,k]; n_gc += dom
            push!(out, @sprintf("    round %d: rank %d node %d  %.3f s%s  gc %.3f s%s  median rank %.4f s", r, k-1, node_of_rank[k], T[r,k],
                  with_solved ? @sprintf(" (%d solves)", Int(m.S[r,k])) : "", G[r,k], dom ? " [mostly GC]" : "", med(T[r,:])))
        end
        nG = count(r -> (k = argmax(T[r,:]); G[r,k] > 0.5*T[r,k]), 1:nr)
        slow = [argmax(T[r,:]) for r in 1:nr]; nodes = unique(node_of_rank)
        push!(out, "    all rounds: slowest rank was >50% GC in $nG of $nr; distinct slowest ranks $(length(unique(slow))) of $np; " *
                   "rounds won by node: " * join(["$nd:$(count(k -> node_of_rank[k]==nd, slow))" for nd in nodes], " "))
    end
    worst(m.C, m.GC1, "compute", true)
    worst(m.X, m.GC2, "exchange", false)
    ratio = [maximum(m.C[r,:]) / max(med(m.C[r,:]), 1e-9) for r in 1:nr]
    push!(out, @sprintf("  slowest/median compute ratio per round: median %.1f, max %.1f; rounds above 10: %d of %d",
                        med(ratio), maximum(ratio), count(>(10), ratio), nr))
    return out
end

function environment_lines(comm, ni)                                      # collective; only rank 0 returns lines
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    envstr = join(["$k=$(get(ENV, k, "<unset>"))" for k in PROFILE_ENV_VARS], ";")
    hv = Float64(hash(envstr) % UInt64(2)^52)
    same_env = MPI.Allreduce(hv, MPI.MIN, comm) == MPI.Allreduce(hv, MPI.MAX, comm)
    cpus_min = MPI.Allreduce(Float64(Sys.CPU_THREADS), MPI.MIN, comm)
    mem_min  = MPI.Allreduce(Sys.total_memory()/1e9, MPI.MIN, comm)
    rank == 0 || return String[]
    counts = [count(==(n), ni.node_of_rank) for n in unique(ni.node_of_rank)]
    gcth = isdefined(Threads, :ngcthreads) ? string(Threads.ngcthreads()) : "n/a"
    mpiv = try strip(first(split(MPI.Get_library_version(), '\n'))) catch; "unavailable" end
    git = try
        c = strip(read(Cmd(`git rev-parse --short HEAD`; dir=@__DIR__), String))
        d = count(!isempty, split(read(Cmd(`git status --porcelain`; dir=@__DIR__), String), '\n'))
        "$c, modified files: $d"
    catch; "unavailable" end
    return ["julia: $VERSION, threads=$(Threads.nthreads()), gc_threads=$gcth, blas_threads=$(LinearAlgebra.BLAS.get_num_threads())",
            "mpi: $mpiv",
            "layout: $np ranks on $(ni.n_nodes) nodes, ranks/node $(minimum(counts))-$(maximum(counts)), min logical CPUs/node $(Int(cpus_min)), min RAM/node $(round(mem_min; digits=0)) GB",
            "env identical on all ranks: $same_env   (as seen by rank 0: " * join(["$k=$(get(ENV, k, "<unset>"))" for k in PROFILE_ENV_VARS], " ") * ")",
            "code: tag=$PROFILE_CODE_TAG src_hash=$PROFILE_SRC_HASH git=$git"]
end

function profile_job(inputs, nparts, distribute, sem, partitioned_model; cfg::RunConfig=config_from_env(), outdir::String=".")
    comm = get_mpi_comm(); rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    #=if get(ENV, "PROFILE_BENCH_PACKED", "false") == "true"
    	bench_packed_cache(comm)
    	return nothing
    end=#
    jobid = get(ENV, "SLURM_JOB_ID", "nojob")
    base = joinpath(outdir, "profile_$(cfg.label)_np$(np)_$(jobid)")
    rep = RunReport(rank == 0 ? open(base * "_report.txt", "w") : nothing)
    ni = node_info(comm)                                                       # collective

    # single-literal @sprintf helpers (the macro needs one string literal as its format)
    f1(x) = @sprintf("%.1f", x);  f2(x) = @sprintf("%.2f", x)
    f3(x) = @sprintf("%.3f", x);  f15(x) = @sprintf("%.15g", x)
    gcint(x) = isfinite(x) ? string(Int(round(x))) : "n/a"

    # ---------------- targets (identical on every rank), then the spatial/angular switches ----------------
    spatial, ang_g, ang_e = generate_stress_test_targets(sem.mesh.gnelem; seed=cfg.seed, frac=cfg.frac, max_e_ext=cfg.max_e_ext)
    cs(v) = Float64(sum(Int(x)*i for (i,x) in enumerate(v); init=0) % 1_000_000_007)
    tv = Float64[length(spatial), cs(spatial), length(ang_g), cs(ang_g), cs(ang_e), sem.mesh.gnelem]
    same_targets = MPI.Allreduce(tv, MPI.MIN, comm) == MPI.Allreduce(tv, MPI.MAX, comm)
    angular_targets = Tuple{Int,Int}[(ang_g[i], ang_e[i]) for i = 1:length(ang_g)]
    cfg.spatial || (spatial = Int[])
    cfg.angular || (angular_targets = Tuple{Int,Int}[])

    # ---------------- setup ----------------
    t0 = time_ns()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=collect(Int, spatial), angular_targets=angular_targets)
    t_setup = MPI.Allreduce((time_ns()-t0)*1e-9, MPI.MAX, comm)
    mesh = setup.mesh

    env_lines = environment_lines(comm, ni)                                    # collective
    axes_flags = let f = [0, 0, 0]; for a in detect_periodic_axes(mesh); f[a] = 1; end; MPI.Allreduce(f, MPI.MAX, comm); end
    _, pk = build_tagged_face_keys_periodic_aware(mesh)
    n_pf = MPI.Allreduce(length(pk), +, comm); n_el = MPI.Allreduce(mesh.nelem, +, comm)
    n_ncf = MPI.Allreduce(mesh.num_ncf, +, comm); n_pairs = MPI.Allreduce(sum(setup.n_dirs), +, comm)
    lv = Float64[mesh.nelem, sum(setup.n_dirs), setup.gid_exch.n_flat, mesh.npoin]
    lmin = MPI.Allreduce(lv, MPI.MIN, comm); lmax = MPI.Allreduce(lv, MPI.MAX, comm); lmean = MPI.Allreduce(lv, MPI.SUM, comm) ./ np

    # ---------------- global switches (set BEFORE the warm-up so warm-up and measured phase run the same code paths) ----------------
    USE_FACE_PLANS[] = cfg.face_plans
    USE_LOCAL_NCF_PLANS[] = cfg.local_plans
    USE_FLAT_FACES[] = cfg.flat_faces
    USE_ANGULAR_PLANS[] = cfg.angular_plans
    USE_UPWIND_ORDER[] = cfg.upwind_order
    USE_INCREMENTAL_PACK[] = cfg.incremental_pack
    USE_GS_INPLACE[] = cfg.gs_inplace
    USE_GS_GC_OFF[] = cfg.gs_gc_off
    GS_MODE[] = cfg.gs_mode
    clear_face_plans!()
    cache_mb = cfg.lu_cache_mb
    budget_lines = String[]
    if cfg.lu_cache_mb < 0                                  # PROFILE_LU_CACHE_MB=-1 means: choose from available memory
    	cache_mb, budget_lines = auto_cache_budget_mb(setup, ni, comm; min_pairs_per_rank=cfg.cache_min_pairs)
    end

    # ---------------- warm-up (always :auto GC) ----------------
    t0 = time_ns()
    if cfg.warmup_outer > 0
        run_unified_manufactured_check_profiled(setup, comm, SweepProfile(); κ_val=cfg.κ_val, σ_val=cfg.σ_val,
            n_outer=cfg.warmup_outer, local_passes=cfg.local_passes, measure_skew=false, coupling=cfg.coupling,
            gc_mode=:auto, gc_every=0, lu_cache_mb=min(cache_mb, 100.0), round_cap = cfg.round_cap, entry_barrier = cfg.entry_barrier)   # small cache: compiles cached AND workspace paths
    end
    t_warm = MPI.Allreduce((time_ns()-t0)*1e-9, MPI.MAX, comm)
    GC.gc(); MPI.Barrier(comm)

    # ---------------- measured phase ----------------
    fill!(ANG_STATS, 0)
    gc0 = gc_snapshot(); prof = SweepProfile(); t0 = time_ns()
    last_delta, I_final, ecache = run_unified_manufactured_check_profiled(setup, comm, prof; κ_val=cfg.κ_val, σ_val=cfg.σ_val,
        n_outer=cfg.n_outer, local_passes=cfg.local_passes, measure_skew=cfg.measure_skew, coupling=cfg.coupling,
        gc_mode=cfg.gc_mode, gc_every=cfg.gc_every, lu_cache_mb=cache_mb, round_cap = cfg.round_cap, entry_barrier = cfg.entry_barrier)
    t_wall = MPI.Allreduce((time_ns()-t0)*1e-9, MPI.MAX, comm)
    gc1 = gc_snapshot()

    red = reduce_profile(prof, comm)                                           # collective
    dgc = gc1 .- gc0
    gc_max = MPI.Allreduce(dgc, MPI.MAX, comm); gc_mean = MPI.Allreduce(dgc, MPI.SUM, comm) ./ np
    rss = MPI.Allreduce(Sys.maxrss()/1e9, MPI.MAX, comm); live = MPI.Allreduce(Base.gc_live_bytes()/1e9, MPI.MAX, comm)
    mats = [gather_sweep_matrices(prof, s, comm) for s in 1:length(prof.sweeps)]      # collective

    # ---------------- every collective diagnostic, same order on all ranks ----------------
    # switch_report MUST come first: verify_local_ncf_plans clears and refills the plan tables, and
    # verify_face_plans builds remote plans, so the counts are only meaningful straight after the measured phase.
    sw_lines     = switch_report(comm)
    gsb_lines    = gs_branch_report(prof, comm)
    diag_lines   = straggler_diagnosis(prof, setup, comm, mats, ni, base)
    lp_lines     = verify_local_ncf_plans(setup, comm)
    ap_lines = verify_angular_plans(setup, I_final, comm)
    verify_lines = verify_face_plans(setup, I_final, comm)
    gt_lines     = verify_gid_table(setup, I_final, comm)
    fp_lines     = solution_fingerprint(mesh, I_final, comm)
    gel_lines    = per_gel_fingerprint_csv(setup, I_final, comm, base * "_gel.csv")
    cache_lines  = elem_cache_report(ecache, sum(setup.n_dirs), comm)

    GC.gc(); MPI.Barrier(comm)                                                 # memory after a full GC (live heap includes the LU cache)
    memv = Float64[Base.gc_live_bytes()/1e9, Sys.maxrss()/1e9, sum(setup.n_dirs), setup.mesh.nelem]
    MEM = reshape(gather_all(memv, comm), 4, np)
    mem_lines = String[]
    if rank == 0
        lo = argmin(MEM[3,:]); hi = argmax(MEM[3,:])
        push!(mem_lines, "memory after a full GC (live heap includes the LU cache): live min $(f3(minimum(MEM[1,:]))) mean $(f3(sum(MEM[1,:])/np)) " *
            "max $(f3(maximum(MEM[1,:]))) GB; peak RSS max $(f3(maximum(MEM[2,:]))) GB; " *
            "fewest pairs: rank $(lo-1) ($(Int(MEM[3,lo])) pairs, $(Int(MEM[4,lo])) elems) live $(f3(MEM[1,lo])) GB; " *
            "most pairs: rank $(hi-1) ($(Int(MEM[3,hi])) pairs, $(Int(MEM[4,hi])) elems) live $(f3(MEM[1,hi])) GB")
    end

    # ---------------- report ----------------
    if rank == 0
        med(v) = sort(v)[(end+1)÷2]
        emit(rep, "##### PROFILE RUN REPORT #####")
        emit(rep, "run: label=$(cfg.label) job=$jobid date=$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        foreach(l -> emit(rep, l), env_lines)
        emit(rep, "config: " * join(["$f=$(getfield(cfg, f))" for f in fieldnames(RunConfig)], " "))
        emit(rep, "gc semantics: warm-up used :auto; measured phase gc_mode=$(cfg.gc_mode) gc_every=$(cfg.gc_every) " *
                  "(off = GC.enable(false) during the measured runner, sync = manual GC.gc every gc_every rounds and between sweeps); " *
                  "GC.gc() + Barrier before the measured phase; no micro-benchmarks ran in this job")
        emit(rep, "problem: gnelem=$(mesh.gnelem) nelem_after=$n_el periodic_axes(x,y,z)=$(Tuple(axes_flags)) periodic_faces=$n_pf " *
                  "same_rank_ncf=$n_ncf (iel,d) pairs=$n_pairs targets(spatial=$(length(spatial)), angular=$(length(angular_targets))) " *
                  "targets identical on all ranks: $same_targets")
        foreach(l -> emit(rep, l), sw_lines)
        emit(rep, "angular mesh: nelemx=$(inputs[:extra_dimensions_nelemx]) nelemy=$(inputs[:extra_dimensions_nelemy]) " *
                  "order=$(inputs[:extra_dimensions_order]) pairs/element=$(f1(n_pairs/mesh.gnelem))")
	effcap(p) = cfg.round_cap < 0 ? max(64, Int(floor(p / (-cfg.round_cap)))) : cfg.round_cap
	emit(rep, "round cap: setting $(cfg.round_cap), effective per-rank cap $(effcap(lmin[2])) to $(effcap(lmax[2])) " *
        	  "pairs per round (0 = unlimited); incremental pack $(cfg.incremental_pack)")
        for (k, name) in enumerate(["elements", "(iel,d) pairs", "exchange entries sent", "CG nodes"])
            emit(rep, "load " * rpad(name, 22) * " min $(Int(lmin[k]))  mean $(f1(lmean[k]))  max $(Int(lmax[k]))  " *
                      "max/mean $(f2(lmax[k]/max(lmean[k], 1e-300)))")
        end

        # correctness block: what must be bit-identical across code changes and rank counts
        foreach(l -> emit(rep, l), fp_lines)
        foreach(l -> emit(rep, l), gel_lines)
        foreach(l -> emit(rep, l), lp_lines)
        foreach(l -> emit(rep, l), ap_lines)
	foreach(l -> emit(rep, l), verify_lines)
        foreach(l -> emit(rep, l), gt_lines)

        # resources and timings
	foreach(l -> emit(rep, l), budget_lines)
        foreach(l -> emit(rep, l), cache_lines)
        foreach(l -> emit(rep, l), mem_lines)
        emit(rep, "timings [s]: setup $(f2(t_setup))  warmup $(f2(t_warm))  measured wall $(f2(t_wall))  " *
                  "sweep_total(max rank) $(f2(red.stat[:sweep_total].tmax))  max_delta after $(cfg.n_outer) outer iterations $(f15(last_delta))")
        emit(rep, "rounds (all measured sweeps): $(Int(red.misc.rounds))   allocated: $(f1(red.misc.alloc/1e9)) GB total   " *
                  "peak RSS(max rank) $(f2(rss)) GB   live heap(max rank) $(f2(live)) GB")
        emit(rep, "gc during measured phase: pause mean $(f3(gc_mean[1])) s max-rank $(f3(gc_max[1])) s; " *
                  "full sweeps mean $(isfinite(gc_mean[2]) ? f1(gc_mean[2]) : "n/a") max $(gcint(gc_max[2]))")

        emit(rep, "section                 mean s     max s       calls")
        for k in (:attempt_solved, :attempt_wasted, :solve, :lu_build, :face_upwind, :classify, :alloc_A, :vol_matrix,
                  :flux_terms, :source_term, :ex_pack, :ex_mpi, :ex_unpack, :allreduce, :skew_exchange, :skew_reduce,
                  :sync_gc, :gs_total, :gs_build, :gs_upwind, :gs_solve, :gs_exchange, :gs_reduce, :empty_pass,
		  :outer_src, :outer_update, :run_fields, :run_init, :run_cache, :periodic_exch, :bookkeeping, :trial_eval, :entry_wait, :entry_wait2)
            s = red.stat[k]
            emit(rep, rpad(String(k), 20) * lpad(f3(s.tsum/np), 11) * lpad(f3(s.tmax), 10) * lpad(string(Int(round(s.nsum))), 12))
        end
        gs_parts = (:gs_build, :gs_upwind, :gs_solve, :gs_exchange, :gs_reduce)
        gs_tot = max(red.stat[:gs_total].tsum / np, 1e-12)
        gs_sum = sum(red.stat[k].tsum for k in gs_parts) / np
        emit(rep, "gs breakdown (mean s per rank, % of gs_total): " *
                  join(["$(replace(String(k), "gs_" => "")) $(f2(red.stat[k].tsum/np)) ($(round(Int, 100*red.stat[k].tsum/np/gs_tot))%)" for k in gs_parts], " | ") *
                  " | untimed $(f2(gs_tot - gs_sum)) ($(round(Int, 100*(gs_tot - gs_sum)/gs_tot))%)")
        foreach(l -> emit(rep, l), gsb_lines)
        for (s, m) in enumerate(mats)
            emit(rep, "sweep $s: $(m.nr) rounds | compute mean-rank $(f2(sum(m.C)/np)) s, sum of round max $(f2(sum(maximum(m.C; dims=2)))) s | " *
                      "exchange mean-rank $(f2(sum(m.X)/np)) s, sum of round max $(f2(sum(maximum(m.X; dims=2)))) s | " *
                      "reduction wait mean-rank $(f2(sum(m.R)/np)) s, max-rank $(f2(maximum(vec(sum(m.R; dims=1))))) s")
        end
        if !isempty(mats)
            emit(rep, "stragglers, sweep 1:")
            foreach(l -> emit(rep, l), straggler_lines(mats[1], ni.node_of_rank))
        end
        foreach(l -> emit(rep, l), diag_lines)

        csvname = base * "_rounds.csv"
        open(csvname, "w") do io
            println(io, "sweep,round,passes_max,solved_sum,compute_max,compute_median,exchange_max,exchange_median,reduce_max,reduce_median,gc_compute_max,gc_exchange_max")
            for (s, m) in enumerate(mats), r in 1:m.nr
                println(io, join(Any[s, r, Int(maximum(m.P[r,:])), Int(sum(m.S[r,:])), maximum(m.C[r,:]), med(m.C[r,:]),
                    maximum(m.X[r,:]), med(m.X[r,:]), maximum(m.R[r,:]), med(m.R[r,:]), maximum(m.GC1[r,:]), maximum(m.GC2[r,:])], ","))
            end
        end
        emit(rep, "files: $(base)_report.txt  $csvname  $(base)_ranktable.csv  $(base)_rankrounds.csv  $(base)_gel.csv")
        emit(rep, "##### END REPORT #####")
        close(rep.io)
    end
    return nothing
end

function build_isotropic_manufactured_fields(setup; I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1, κ_val=0.001, σ_val=0.0005)
    mesh = setup.mesh; ngl = mesh.ngl
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)
    κ = fill(κ_val, mesh.npoin)

    groups, representative_of = group_elements_by_angular_mesh(mesh)
    coeff_all = Dict{Int,Vector{Float64}}(); Ωm_all = Dict{Int,Vector{NTuple{3,Float64}}}()
    for rep in keys(groups)
        em = mesh.extra_mesh[rep]
        w_ext = zeros(em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                w_ext[em.extra_connijk[e_ext,iθ,iϕ]] += em.extra_metrics.Je[e_ext,iθ,iϕ]*em.ωθ[iθ]*em.ωθ[iϕ]
            end
        end
        coeff_all[rep] = σ_val .* w_ext ./ sum(w_ext)
        Ωm_all[rep] = [(sin(em.extra_coords[1,m])*cos(em.extra_coords[2,m]),
                        sin(em.extra_coords[1,m])*sin(em.extra_coords[2,m]),
                        cos(em.extra_coords[1,m])) for m = 1:em.extra_npoin]
    end

    q_ext_all = Dict{Tuple{Int,Int},Vector{Float64}}()          # LOCAL format: ngl^3 entries in lidx order
    bc_fn_all = Dict{Tuple{Int,Int},Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        rep = representative_of[iel]; coeff = coeff_all[rep]; Ω_m = Ωm_all[rep]
        scat = zeros(ngl^3)                                      # direction-independent for isotropic scattering
        for k=1:ngl, j=1:ngl, i=1:ngl
            gip = mesh.connijk[iel,i,j,k]
            s = 0.0
            for m = 1:length(coeff)
                s += coeff[m] * I_manuf(mesh.x[gip], mesh.y[gip], mesh.z[gip], Ω_m[m]...)
            end
            scat[lidx(i,j,k,ngl)] = s
        end
        for dk = 1:em.extra_npoin
            θk = em.extra_coords[1,dk]; ϕk = em.extra_coords[2,dk]
            Ωxk = sin(θk)*cos(ϕk); Ωyk = sin(θk)*sin(ϕk); Ωzk = cos(θk)
            q_vec = zeros(ngl^3)
            for k=1:ngl, j=1:ngl, i=1:ngl
                gip = mesh.connijk[iel,i,j,k]; p = lidx(i,j,k,ngl)
                x_,y_,z_ = mesh.x[gip], mesh.y[gip], mesh.z[gip]
                grad = dI_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                q_vec[p] = Ωxk*grad[1]+Ωyk*grad[2]+Ωzk*grad[3] + κ_val*I_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk) - scat[p]
            end
            q_ext_all[(iel,dk)] = q_vec
            bc_fn_all[(iel,dk)] = (x,y,z,ip,nx,ny,nz) -> I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
        end
    end
    return (κ=κ, bc_fn_all=bc_fn_all, q_ext_all=q_ext_all, coeff_all=coeff_all, representative_of=representative_of)
end

function assemble_isotropic_source(mesh, coeff_all, representative_of, q_ext_all, I_cg_all, variant::Symbol)
    ngl = mesh.ngl
    out = Dict{Tuple{Int,Int},Vector{Float64}}()
    for iel = 1:mesh.nelem
        nd = mesh.extra_mesh[iel].extra_npoin
        coeff = coeff_all[representative_of[iel]]
        if variant == :legacy
            for dk = 1:nd
                q = copy(q_ext_all[(iel,dk)])
                for m = 1:nd
                    coeff[m] == 0.0 && continue
                    q .+= coeff[m] .* I_cg_all[(iel,m)]
                end
                out[(iel,dk)] = q
            end
        else
            Imat = Matrix{Float64}(undef, ngl^3, nd)
            for m = 1:nd; Imat[:,m] .= I_cg_all[(iel,m)]; end
            if variant == :isotropic
                gsrc = Imat * coeff
                for dk = 1:nd; out[(iel,dk)] = q_ext_all[(iel,dk)] .+ gsrc; end
            elseif variant == :matprod
                C = repeat(reshape(coeff, 1, :), nd, 1)          # C[dk,m] = coeff[m]; general phase functions fill C fully
                S = Imat * transpose(C)
                for dk = 1:nd; out[(iel,dk)] = q_ext_all[(iel,dk)] .+ view(S, :, dk); end
            else
                error("unknown coupling variant $variant")
            end
        end
    end
    return out
end

function run_unified_manufactured_check_profiled(setup, comm, prof::SweepProfile;
        κ_val=0.001, σ_val=0.0005, n_outer::Int=2, local_passes::Bool=true,
        measure_skew::Bool=true, coupling::Symbol=:isotropic, lu_cache_mb::Float64=0.0, gc_mode::Symbol=:auto, gc_every::Int=8,
        round_cap::Int=0, entry_barrier::Bool=false, field_kw...)
    mesh = setup.mesh; metrics = setup.metrics
    ψ, dψ, ω = setup.sem.basis.ψ, setup.sem.basis.dψ, setup.sem.ω
    ngl = mesh.ngl; n_dirs = setup.n_dirs
    t0 = ptic()
    f = build_isotropic_manufactured_fields(setup; κ_val=κ_val, σ_val=σ_val, field_kw...)
    ptoc!(prof, :run_fields, t0)

    t0 = ptic()
    I_cg_all = Dict{Tuple{Int,Int},Vector{Float64}}()
    for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]; I_cg_all[(iel,dd)] = zeros(ngl^3); end
    ptoc!(prof, :run_init, t0)

    last_delta = NaN

    gc_prev = gc_mode == :auto ? true : GC.enable(false)
    t0 = ptic()
    cache = make_elem_cache(mesh, metrics, n_dirs, ngl; budget_mb=lu_cache_mb)
    ptoc!(prof, :run_cache, t0)

    try
    	for it = 1:n_outer
            t0 = ptic()
            extra_source_all = assemble_isotropic_source(mesh, f.coeff_all, f.representative_of, f.q_ext_all, I_cg_all, coupling)
            ptoc!(prof, :outer_src, t0)
            if entry_barrier
    		tb = time_ns(); MPI.Barrier(comm)
    		padd!(prof, it == 1 ? :entry_wait : :entry_wait2, (time_ns() - tb) * 1e-9)
	    end
	    x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            	setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence, setup.periodic_face_correspondence,
            	setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group, setup.neighbor_is_periodic,
            	setup.same_rank_kind, setup.same_rank_partner, setup.same_rank_is_periodic,
            	setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
            	setup.remote_gid_coords,
            	setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh,
            	setup.gid_exch, setup.periodic_gid_exch,
            	f.κ, f.bc_fn_all, comm;
            	prof=prof, local_passes=local_passes, elem_cache=cache, measure_skew=measure_skew,
            	extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all, gc_mode=gc_mode, gc_every=gc_every, round_cap=round_cap)

            t0 = ptic()
            dl = 0.0
            for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
            	dl = max(dl, maximum(abs.(x_dg[(iel,dd)] .- I_cg_all[(iel,dd)])))
            	I_cg_all[(iel,dd)] = copy(x_dg[(iel,dd)])
            end
            last_delta = MPI.Allreduce(dl, MPI.MAX, comm)
            ptoc!(prof, :outer_update, t0)
            if gc_mode == :sync                                  # also collect between sweeps
            	GC.enable(true); GC.gc(false); GC.enable(false)
            end
	end
    finally
    	gc_mode == :auto || GC.enable(gc_prev)
    end
    return last_delta, I_cg_all, cache
end

function bench_allreduce(comm; n::Int=200)
    x = Int32[1,2]
    MPI.Allreduce(x, +, comm); MPI.Barrier(comm)                           # warm-up
    t0 = time_ns(); for _ = 1:n; MPI.Allreduce(x, +, comm); end
    t_ar = MPI.Allreduce((time_ns()-t0)*1e-9/n, MPI.MAX, comm)
    MPI.Barrier(comm)
    t0 = time_ns(); for _ = 1:n; MPI.Barrier(comm); end
    t_b = MPI.Allreduce((time_ns()-t0)*1e-9/n, MPI.MAX, comm)
    MPI.Comm_rank(comm) == 0 && println("bare collectives, $(MPI.Comm_size(comm)) ranks: " *
        "Allreduce(2 ints) = $(round(1e6*t_ar; digits=1)) us, Barrier = $(round(1e6*t_b; digits=1)) us")
end

function bench_coupling_variants(setup, comm; reps::Int=2, field_kw...)
    mesh = setup.mesh; ngl = mesh.ngl
    f = build_isotropic_manufactured_fields(setup; field_kw...)
    I = Dict{Tuple{Int,Int},Vector{Float64}}(
        (iel,d) => [1.0 + 1e-3*(iel + d + p) for p = 1:ngl^3] for iel = 1:mesh.nelem for d = 1:setup.n_dirs[iel])
    res = Dict{Symbol,Any}(); tm = Dict{Symbol,Float64}()
    for v in (:isotropic, :matprod, :legacy)
        res[v] = assemble_isotropic_source(mesh, f.coeff_all, f.representative_of, f.q_ext_all, I, v)       # warm-up
        t = Inf
        for _ = 1:reps
            t0 = time_ns()
            res[v] = assemble_isotropic_source(mesh, f.coeff_all, f.representative_of, f.q_ext_all, I, v)
            t = min(t, (time_ns()-t0)*1e-9)
        end
        tm[v] = MPI.Allreduce(t, MPI.MAX, comm)
    end
    dmax(a, b) = (m = 0.0; for (k,v) in a; m = max(m, maximum(abs.(v .- b[k]))); end; m)
    d_mp = MPI.Allreduce(dmax(res[:matprod], res[:isotropic]), MPI.MAX, comm)
    d_lg = MPI.Allreduce(dmax(res[:legacy],  res[:isotropic]), MPI.MAX, comm)
    MPI.Comm_rank(comm) == 0 && println("coupling assembly, max rank [s]: isotropic=$(round(tm[:isotropic]; sigdigits=3)) " *
        "matprod=$(round(tm[:matprod]; sigdigits=3)) legacy=$(round(tm[:legacy]; sigdigits=3));  " *
        "max diff vs isotropic: matprod=$d_mp legacy=$d_lg")
end

function bench_local_solves(setup, comm; κ_val::Float64=0.001, n_samples::Int=24, reps::Int=40)
    rank = MPI.Comm_rank(comm)
    mesh = setup.mesh; metrics = setup.metrics
    ψ, dψ, ω = setup.sem.basis.ψ, setup.sem.basis.dψ, setup.sem.ω
    ngl = mesh.ngl; n3 = ngl^3
    κ = fill(κ_val, mesh.npoin)
    mats = Matrix{Float64}[]
    for s = 1:n_samples
        iel = 1 + mod(37*s, mesh.nelem)
        d = 1 + mod(53*s, mesh.extra_mesh[iel].extra_npoin)
        em = mesh.extra_mesh[iel]; θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz)
        A = zeros(n3,n3); rhs = zeros(n3)
        build_dg_volume_matrix!(A, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
            metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz,
            metrics.Je, κ, mesh.connijk)
        for fid = 1:6
            is_active[fid] || continue
            fn = compute_face_normal_avg(iel, fid, ngl,
                metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz)
            build_dg_flux_terms!(A, rhs, fid, ngl, Ωx,Ωy,Ωz, is_inflow[fid], is_active[fid], fn, ω,
                metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz,
                metrics.Je, iel, (i,j,k)->0.0)
        end
        push!(mats, A)
    end
    b = [1.0 + 0.01*p for p = 1:n3]
    function timeit(f)
        f()                                                       # warm-up
        t0 = time_ns(); for _ = 1:reps; f(); end
        return (time_ns()-t0)*1e-9/reps
    end
    t_bs = 0.0; t_lu = 0.0; t_ld = 0.0; t_slu = 0.0; t_sld = 0.0; nnzA = 0.0; nnzF = 0.0; xres = 0.0
    x = similar(b)
    for A in mats
        t_bs += timeit(() -> A \ b)
        t_lu += timeit(() -> lu(A))
        F = lu(A)
        t_ld += timeit(() -> (x .= b; ldiv!(F, x)))
        xres = max(xres, maximum(abs.(A*x .- b)))
        S = sparse(A); nnzA += nnz(S)
        t_slu += timeit(() -> lu(S))
        SF = lu(S); nnzF += nnz(SF.L) + nnz(SF.U)
        t_sld += timeit(() -> SF \ b)
    end
    n = length(mats); us(t) = round(1e6*t/n; digits=1)
    pairs_global = MPI.Allreduce(sum(setup.n_dirs), +, comm)
    if rank == 0
        println("local solves, $n real element matrices, $(n3)x$(n3):")
        println("  matrix nnz: $(round(nnzA/n; digits=0)) ($(round(100*nnzA/n/n3^2; digits=1))% dense);  residual check max|Ax-b| = $xres")
        println("  dense   A\\b = $(us(t_bs)) us | lu(A) = $(us(t_lu)) us | ldiv!(F,b) = $(us(t_ld)) us  -> cached-LU speedup on the solve alone: $(round(t_bs/t_ld; digits=1))x")
        println("  sparse  lu = $(us(t_slu)) us | solve = $(us(t_sld)) us | factor nnz = $(round(nnzF/n; digits=0))")
        println("  cache memory for $pairs_global pairs: dense = $(round(pairs_global*n3^2*8/1e9; digits=1)) GB, " *
                "sparse LU ~ $(round(pairs_global*nnzF/n*12/1e9; digits=1)) GB (all ranks together)")
    end
end

function reduce_profile(prof::SweepProfile, comm)
    np = MPI.Comm_size(comm)
    stat = Dict{Symbol,Any}()
    for k in PROF_SECTIONS
        t = get(prof.t, k, 0.0); n = Float64(get(prof.n, k, 0))
        stat[k] = (tmax=MPI.Allreduce(t, MPI.MAX, comm), tmin=MPI.Allreduce(t, MPI.MIN, comm),
                   tsum=MPI.Allreduce(t, MPI.SUM, comm), nsum=MPI.Allreduce(n, MPI.SUM, comm))
    end
    misc = (gc=MPI.Allreduce(get(prof.t, :gc_seconds, 0.0), MPI.MAX, comm),
            alloc=MPI.Allreduce(Float64(get(prof.n, :alloc_bytes, 0)), MPI.SUM, comm),
            rounds=MPI.Allreduce(Float64(get(prof.n, :rounds, 0)), MPI.MAX, comm), np=np)
    red(v, op) = MPI.Allreduce(Float64.(v), op, comm)
    sweeps = NamedTuple[]
    for L in prof.sweeps
        push!(sweeps, (nr=length(L.passes), solved=red(L.solved, +), passes_max=red(L.passes, MPI.MAX),
                       t_compute=red(L.t_compute, MPI.MAX), t_exchange=red(L.t_exchange, MPI.MAX),
                       t_skew=red(L.t_skew, MPI.MAX), t_reduce=red(L.t_reduce, MPI.MAX),
                       sent=red(L.sent, +), valid=red(L.valid, +), fresh=red(L.fresh, +), t_empty=red(L.t_empty, MPI.MAX)))
    end
    return (stat=stat, misc=misc, sweeps=sweeps)
end

function print_sweep_profile(red, comm, ngl::Int; label::String, n_outer::Int, t_total::Float64, last_delta::Float64)
    MPI.Comm_rank(comm) == 0 || return
    stat, misc, sweeps = red.stat, red.misc, red.sweeps
    np = misc.np; ref = max(stat[:sweep_total].tmax, 1e-12)
    println("\n================ SWEEP PROFILE [$label] ($np ranks, $n_outer outer iterations, ngl=$ngl) ================")
    println("BLAS threads/rank: $(LinearAlgebra.BLAS.get_num_threads()), Julia threads: $(Threads.nthreads())")
    println("runner wall time (max rank): $(round(t_total; digits=2)) s;  sweep_total (max rank): $(round(ref; digits=2)) s;  " *
            "max_delta after $n_outer outer iterations: $last_delta (not converged on purpose)")
    println("communication rounds (all sweeps): $(Int(misc.rounds));  GC: $(round(misc.gc; digits=2)) s max;  allocated: $(round(misc.alloc/1e9; digits=2)) GB total")
    println(rpad("section", 16), lpad("max s", 10), lpad("min s", 10), lpad("mean s", 10), lpad("max/mean", 10),
            lpad("calls", 12), lpad("us/call", 10), lpad("% sweep", 9))
    for k in PROF_SECTIONS
        s = stat[k]; mean_t = s.tsum/np
        println(rpad(String(k), 16), lpad(round(s.tmax; digits=3), 10), lpad(round(s.tmin; digits=3), 10),
                lpad(round(mean_t; digits=3), 10), lpad(round(s.tmax/max(mean_t,1e-12); digits=2), 10),
                lpad(Int(round(s.nsum)), 12), lpad(round(1e6*s.tsum/max(s.nsum,1.0); digits=1), 10),
                lpad(round(100*s.tmax/ref; digits=1), 9))
    end
    println("(attempt_solved/attempt_wasted contain the component rows; empty_pass is a subset of attempt_wasted; " *
            "outer_src/outer_update are outside sweep_total)")
    ex_mean = (stat[:ex_pack].tsum + stat[:ex_mpi].tsum + stat[:ex_unpack].tsum)/np
    println("exchange, mean over ranks: pack $(round(stat[:ex_pack].tsum/np; digits=3)) s + mpi $(round(stat[:ex_mpi].tsum/np; digits=3)) s + " *
            "unpack $(round(stat[:ex_unpack].tsum/np; digits=3)) s = $(round(ex_mean; digits=3)) s")
    att_s = stat[:attempt_solved].nsum; att_w = stat[:attempt_wasted].nsum
    println("attempts: $(Int(att_s)) solved, $(Int(att_w)) blocked ($(round(100*att_w/max(att_s+att_w,1); digits=1))% of attempts), " *
            "blocked cost $(round(100*stat[:attempt_wasted].tsum/max(stat[:attempt_solved].tsum+stat[:attempt_wasted].tsum,1e-12); digits=1))% of attempt time; " *
            "empty passes: $(Int(stat[:empty_pass].nsum)) calls, $(round(stat[:empty_pass].tsum/np; digits=3)) s mean per rank")
    n3 = ngl^3; ns = stat[:solve].nsum
    if ns > 0 && stat[:solve].tsum > 0
        println("dense solves: $(Int(ns)) of size $n3 -> $(round(1e6*stat[:solve].tsum/ns; digits=1)) us each, ~" *
                "$(round(((2/3)*n3^3 + 2*n3^2)*ns/stat[:solve].tsum/1e9; digits=2)) GFLOP/s per rank")
    end

    tot_sent = 0.0; tot_fresh = 0.0
    for (s, L) in enumerate(sweeps)
        tot_sent += sum(L.sent); tot_fresh += sum(L.fresh)
        println("  sweep $s: $(L.nr) rounds, passes/round (max over ranks) mean $(round(sum(L.passes_max)/L.nr; digits=2)) max $(Int(maximum(L.passes_max))), " *
                "time (max rank) compute $(round(sum(L.t_compute); digits=2)) s, exchange $(round(sum(L.t_exchange); digits=2)) s, " *
                "skew waits $(round(sum(L.t_skew); digits=2)) s, reduction $(round(sum(L.t_reduce); digits=2)) s")
    end
    println("exchange redundancy over all rounds: $(round(100*tot_fresh/max(tot_sent,1); digits=1))% of sent entries were new " *
            "(a changed-only exchange would send about that fraction)")
    if !isempty(sweeps)
        L = sweeps[1]; tot = max(sum(L.solved), 1.0); cum = 0.0
        println("\nper-round table, sweep 1 (solved/sent are summed over ranks; times are the max over ranks):")
        println(lpad("round",6), lpad("solved",10), lpad("cum%",7), lpad("passes",8), lpad("compute s",11), lpad("exchange s",11),
                lpad("skew s",9), lpad("reduce s",10), lpad("sent",10), lpad("valid%",8), lpad("new%",7))
        for r = 1:L.nr
            cum += L.solved[r]
            println(lpad(r,6), lpad(Int(L.solved[r]),10), lpad(round(100*cum/tot; digits=1),7), lpad(Int(L.passes_max[r]),8),
                    lpad(round(L.t_compute[r]; digits=3),11), lpad(round(L.t_exchange[r]; digits=3),11),
                    lpad(round(L.t_skew[r]; digits=3),9), lpad(round(L.t_reduce[r]; digits=3),10),
                    lpad(Int(L.sent[r]),10), lpad(round(100*L.valid[r]/max(L.sent[r],1); digits=1),8),
                    lpad(round(100*L.fresh[r]/max(L.sent[r],1); digits=1),7))
        end
    end
end

function print_mode_comparison(a, b, comm)
    MPI.Comm_rank(comm) == 0 || return
    np = a.misc.np
    println("\n================ single_pass vs local_passes ================")
    println(rpad("quantity", 44), lpad("single_pass", 14), lpad("local_passes", 14))
    row(name, f) = println(rpad(name, 44), lpad(round(f(a); sigdigits=4), 14), lpad(round(f(b); sigdigits=4), 14))
    row("communication rounds", r -> r.misc.rounds)
    row("sweep_total, max rank [s]", r -> r.stat[:sweep_total].tmax)
    row("solved attempts, mean [s]", r -> r.stat[:attempt_solved].tsum/np)
    row("blocked attempts [calls]", r -> r.stat[:attempt_wasted].nsum)
    row("blocked attempts, mean [s]", r -> r.stat[:attempt_wasted].tsum/np)
    row("empty passes [calls]", r -> r.stat[:empty_pass].nsum)
    row("exchange (pack+mpi+unpack), mean [s]", r -> (r.stat[:ex_pack].tsum + r.stat[:ex_mpi].tsum + r.stat[:ex_unpack].tsum)/np)
    row("allreduce, max rank [s]", r -> r.stat[:allreduce].tmax)
    row("skew waits, max rank [s]", r -> r.stat[:skew_exchange].tmax + r.stat[:skew_reduce].tmax)
end

function test_profile_stress_unified_manufactured(inputs, nparts, distribute, sem, partitioned_model;
        seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4, n_outer::Int=2, warmup::Bool=true,
        compare_modes::Bool=true, coupling::Symbol=:isotropic, measure_skew::Bool=true,
        run_micro_benchmarks::Bool=true, κ_val=0.001, σ_val=0.0005)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()

    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(sem.mesh.gnelem; seed=seed, frac=frac, max_e_ext=max_e_ext)
    check_targets_consistent(spatial_target_gels, angular_target_gels, angular_target_e_exts, sem.mesh.gnelem, comm)
    angular_targets = Tuple{Int,Int}[(angular_target_gels[i], angular_target_e_exts[i]) for i = 1:length(angular_target_gels)]

    bench_allreduce(comm)
    bench_jitter(comm)

    t0 = time_ns()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=collect(Int, spatial_target_gels), angular_targets=angular_targets)
    t_setup = MPI.Allreduce((time_ns() - t0)*1e-9, MPI.MAX, comm)
    report_rank_load(setup, comm)
    profile_gc_modes(setup, comm) 
    report_neighbors(setup, comm, [42, 43, 9, 44, 45, 62, 17])

    rank == 0 && println("setup (refinement + all tables + exchanges): $(round(t_setup; digits=2)) s")
    print_run_fingerprint(setup.mesh, length(spatial_target_gels), length(angular_targets), comm; label="profile")

    if run_micro_benchmarks
        bench_coupling_variants(setup, comm; κ_val=κ_val, σ_val=σ_val)
        bench_local_solves(setup, comm; κ_val=κ_val)
    end

    if warmup                                              # compile everything once so JIT time doesn't land in the profile
        rank == 0 && println("warm-up iteration ...")
        run_unified_manufactured_check_profiled(setup, comm, SweepProfile(); κ_val=κ_val, σ_val=σ_val,
            n_outer=1, local_passes=true, measure_skew=false, coupling=coupling)
    end

    modes = compare_modes ? [(:single_pass, false), (:local_passes, true)] : [(:local_passes, true)]
    reduced = Dict{Symbol,Any}(); finals = Dict{Symbol,Any}()
    for (label, lp) in modes
        rank == 0 && println("\n##### profiling mode: $label #####")
        prof = SweepProfile()
        t0 = time_ns()
        last_delta, I_final = run_unified_manufactured_check_profiled(setup, comm, prof; κ_val=κ_val, σ_val=σ_val,
            n_outer=n_outer, local_passes=lp, measure_skew=measure_skew, coupling=coupling)
        t_total = MPI.Allreduce((time_ns() - t0)*1e-9, MPI.MAX, comm)
        red = reduce_profile(prof, comm)
        print_sweep_profile(red, comm, setup.mesh.ngl; label=String(label), n_outer=n_outer, t_total=t_total, last_delta=last_delta)
        report_stragglers(prof, comm)
	reduced[label] = red; finals[label] = I_final
    end

    if compare_modes
        dmax = 0.0
        for (k, v) in finals[:single_pass]
            dmax = max(dmax, maximum(abs.(v .- finals[:local_passes][k])))
        end
        dmax = MPI.Allreduce(dmax, MPI.MAX, comm)
        rank == 0 && println("\nmax |I(single_pass) - I(local_passes)| = $dmax  " *
                              (dmax == 0.0 ? "✓ bit-identical" : dmax < 1e-12 ? "✓ roundoff" : "✗ the loop change altered the answer"))
        print_mode_comparison(reduced[:single_pass], reduced[:local_passes], comm)
    end
    return reduced
end

function print_run_fingerprint(mesh, n_spatial_targets, n_angular_targets, comm; label="run")
    rank = MPI.Comm_rank(comm)
    flags = [0, 0, 0]
    for a in detect_periodic_axes(mesh); flags[a] = 1; end
    flags = MPI.Allreduce(flags, MPI.MAX, comm)
    _, pk = build_tagged_face_keys_periodic_aware(mesh)
    n_pf  = MPI.Allreduce(length(pk), +, comm)
    n_el  = MPI.Allreduce(mesh.nelem, +, comm)
    n_ncf = MPI.Allreduce(mesh.num_ncf, +, comm)         # same-rank interfaces only, so it varies with the partition
    rank == 0 && println("[fingerprint $label] ranks=$(MPI.Comm_size(comm)) gnelem=$(mesh.gnelem) " *
        "nelem_after=$n_el periodic_axes(x,y,z)=$(Tuple(flags)) periodic_faces=$n_pf " *
        "targets(spatial=$n_spatial_targets, angular=$n_angular_targets) same_rank_ncf=$n_ncf")
end

function check_refined_set_periodic_closure(mesh, ncf_as_child, ncf_as_parent, gel_to_class, comm)
    rank = MPI.Comm_rank(comm); ngl = mesh.ngl
    function diam(iel)
        lo = fill(Inf, 3); hi = fill(-Inf, 3)
        for k=1:ngl, j=1:ngl, i=1:ngl
            ip = mesh.connijk[iel,i,j,k]; c = (mesh.x[ip], mesh.y[ip], mesh.z[ip])
            for a = 1:3; lo[a] = min(lo[a], c[a]); hi[a] = max(hi[a], c[a]); end
        end
        return maximum(hi .- lo)
    end
    sizes = [diam(i) for i = 1:mesh.nelem]
    base_h = MPI.Allreduce(maximum(sizes; init=0.0), MPI.MAX, comm)
    gather(v::Vector{Float64}) = (cnt = MPI.Allgather(Int32(length(v)), comm);
                                  Set{Int}(Int(round(x)) for x in MPI.Allgatherv(v, Cint.(cnt), comm)))
    fine = gather(Float64[Float64(mesh.el2gel[i]) for i = 1:mesh.nelem if sizes[i] < 0.75*base_h])
    tl = Float64[]
    for (iel,_) in keys(ncf_as_child);  push!(tl, Float64(mesh.el2gel[iel])); end
    for (iel,_) in keys(ncf_as_parent); push!(tl, Float64(mesh.el2gel[iel])); end
    touching = gather(tl)
    gap = 0; n_buffer = 0
    for g in fine, c in get(gel_to_class, g, ())
        Int(c) in fine || (gap += 1)                         # fine element with a coarse class partner: a real hole
    end
    for g in touching, c in get(gel_to_class, g, ())
        (g in fine) || (Int(c) in touching) || (n_buffer += 1)   # coarse buffer element, partner not near refinement: expected
    end
    rank == 0 && println("refined-element periodic closure: ", gap == 0 ? "✓" : "✗ $gap fine elements have a coarse class partner",
                          ";  expected coarse-buffer asymmetries: $n_buffer")
    return gap == 0
end

function local_gips_lidx(mesh, iel)
    ngl = mesh.ngl
    g = Vector{eltype(mesh.connijk)}(undef, ngl^3)
    for k=1:ngl, j=1:ngl, i=1:ngl
        g[lidx(i,j,k,ngl)] = mesh.connijk[iel,i,j,k]
    end
    return g
end

function check_targets_consistent(spatial_gels, ang_gels, ang_exts, gnelem, comm)
    rank = MPI.Comm_rank(comm)
    cs(v) = sum(Int(x)*i for (i,x) in enumerate(v); init=0)       # order-sensitive checksum
    vals = Int[length(spatial_gels), cs(spatial_gels), length(ang_gels), cs(ang_gels), cs(ang_exts), Int(gnelem)]
    same = MPI.Allreduce(vals, MPI.MIN, comm) == MPI.Allreduce(vals, MPI.MAX, comm)
    rank == 0 && println("targets identical on all ranks: $same  (n_spatial=$(vals[1]), n_angular=$(vals[3]), gnelem=$(vals[6]))")
    return same
end

function gather_all(v::Vector{Float64}, comm)
    counts = MPI.Allgather(Int32(length(v)), comm)
    return MPI.Allgatherv(v, Cint.(counts), comm)            # concatenated in rank order
end

function node_info(comm)
    rank = MPI.Comm_rank(comm)
    host = gethostname()
    h = 0
    for (i, c) in enumerate(host); h = mod(h*131 + Int(c) + i, 1_073_741_789); end    # deterministic, no salting
    hs = Int.(round.(gather_all(Float64[h], comm)))
    uniq = sort(unique(hs))
    return (node_id = findfirst(==(h), uniq),
            local_rank = count(==(h), hs[1:rank]),            # ranks below me on the same node
            ranks_on_node = count(==(h), hs),
            node_of_rank = [findfirst(==(x), uniq) for x in hs],
            n_nodes = length(uniq), host = host)
end

function report_rank_load(setup, comm)
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    v = Float64[setup.mesh.nelem, sum(setup.n_dirs), setup.gid_exch.n_flat, setup.mesh.npoin]
    mn = MPI.Allreduce(v, MPI.MIN, comm); mx = MPI.Allreduce(v, MPI.MAX, comm); sm = MPI.Allreduce(v, MPI.SUM, comm)
    if rank == 0
        for (k, name) in enumerate(["elements", "(iel,d) pairs", "exchange entries sent", "CG nodes"])
            println(rpad(name, 24), "min ", Int(mn[k]), "  mean ", round(sm[k]/np; digits=1), "  max ", Int(mx[k]),
                    "  max/mean ", round(mx[k]/(sm[k]/np); digits=2))
        end
    end
end

function bench_jitter(comm; seconds::Float64=12.0, n3::Int=125, chunk_solves::Int=50)
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    ni = node_info(comm)
    A = rand(n3, n3) + n3*I; b = rand(n3); x = similar(b)
    chunk() = for _ = 1:chunk_solves; x .= A \ b; end
    chunk(); chunk()

    function phase(active::Bool)
        MPI.Barrier(comm)
        ts = Float64[]; gcs = Float64[]
        if active
            t_end = time_ns() + round(UInt64, seconds*1e9)
            while time_ns() < t_end
                g0 = Base.gc_time_ns(); t0 = time_ns()
                chunk()
                push!(ts, (time_ns() - t0)*1e-9); push!(gcs, (Base.gc_time_ns() - g0)*1e-9)
            end
        else
            sleep(seconds)
        end
        MPI.Barrier(comm)
        return ts, gcs
    end
    function stats(ts, gcs)                        # median, worst, spikes (>3x own median), chunks, GC-dominated spikes
        isempty(ts) && return Float64[NaN, NaN, 0.0, 0.0, 0.0]
        m = sort(ts)[(length(ts)+1)÷2]
        sp = [i for i in eachindex(ts) if ts[i] > 3*m]
        return Float64[m, maximum(ts), length(sp), length(ts), count(i -> gcs[i] > 0.5*ts[i], sp)]
    end

    tsA, gcA = phase(true)
    tsB, gcB = phase(ni.local_rank == 0)
    gA = gather_all(stats(tsA, gcA), comm); gB = gather_all(stats(tsB, gcB), comm)
    if rank == 0
        SA = reshape(gA, 5, np); SB = reshape(gB, 5, np); N = ni.node_of_rank
        counts = [count(==(n), N) for n in unique(N)]
        us(t) = round(1e6*t/chunk_solves; digits=1)
        medsA = sort(SA[1,:]); medA = medsA[(end+1)÷2]
        medsB = sort(filter(!isnan, SB[1,:])); medB = medsB[(end+1)÷2]
        println("node layout: $(ni.n_nodes) nodes, ranks/node min $(minimum(counts)) max $(maximum(counts)); " *
                "logical CPUs on rank 0's node: $(Sys.CPU_THREADS)")
        println("phase A (all ranks): per-solve median over ranks $(us(medA)) us, slowest rank's median $(us(medsA[end])) us; " *
                "spikes: $(Int(sum(SA[3,:]))) in $(count(>(0), SA[3,:])) ranks, of which GC-dominated: $(Int(sum(SA[5,:]))); " *
                "worst chunk $(round(1e3*maximum(SA[2,:]); digits=1)) ms (typical $(round(1e3*medA; digits=1)) ms)")
        println("phase B (one rank per node): per-solve median $(us(medB)) us  ->  contention factor A/B = $(round(medA/medB; digits=2))")
        for r in sortperm(SA[2,:]; rev=true)[1:min(5,np)]
            println("  worst chunk: rank $(r-1) (node $(N[r])): $(round(1e3*SA[2,r]; digits=1)) ms, own median $(round(1e3*SA[1,r]; digits=1)) ms")
        end
    end
end

function solution_fingerprint(mesh, I_cg_all, comm)
    s = 0.0; c = 0.0; mx = 0.0; n = 0
    for iel = 1:mesh.nelem
        gel = Float64(mesh.el2gel[iel]); em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            v = I_cg_all[(iel,d)]
            for p in eachindex(v)
                y = sin(0.7*gel + 1.3*d + 0.017*p) * v[p] - c      # Kahan-compensated weighted sum
                t = s + y; c = (t - s) - y; s = t
                mx = max(mx, abs(v[p])); n += 1
            end
        end
    end
    tot = MPI.Allreduce(Float64[s, n], +, comm); m = MPI.Allreduce(mx, MPI.MAX, comm)
    return MPI.Comm_rank(comm) == 0 ?
        String[@sprintf("solution fingerprint: n=%d weighted_sum=%.15g max=%.15g", Int(tot[2]), tot[1], m)] : String[]
end


function gs_branch_report(prof::SweepProfile, comm)
    np = MPI.Comm_size(comm); nb = length(BRANCH_LABELS)
    v = Float64[get(prof.t, :gs_upwind, 0.0)]
    for b in BRANCH_LABELS; push!(v, get(prof.t, GS_BRANCH_UP[b], 0.0)); end
    for b in BRANCH_LABELS; push!(v, get(prof.t, GS_BRANCH_EV[b], 0.0)); end
    for b in BRANCH_LABELS; push!(v, Float64(get(prof.n, GS_BRANCH_UP[b], 0))); end
    M = reshape(gather_all(v, comm), length(v), np)
    MPI.Comm_rank(comm) == 0 || return String[]
    w = argmax(M[1,:])
    out = String["GS face work by branch on rank $(w-1), the rank with the largest gs_upwind " *
                 "($(round(M[1,w]; digits=2)) s; mean over ranks $(round(sum(M[1,:])/np; digits=2)) s):"]
    for (j, b) in enumerate(BRANCH_LABELS)
        calls = M[1+2nb+j, w]; calls == 0 && continue
        up = M[1+j, w]; ev = M[1+nb+j, w]
        push!(out, "  $(rpad(String(b), 16)) construct $(round(up; digits=3)) s ($(round(1e6*up/calls; digits=1)) us/call), " *
                   "evaluate $(round(ev; digits=3)) s, calls $(Int(calls))")
    end
    return out
end

function switch_report(comm)      # REPLACES the existing definition
    n = Float64[sum(length(d) for d in values(LOCAL_FINE_PLANS); init=0),
                sum(length(d) for d in values(LOCAL_COARSE_PLANS); init=0),
                sum(length(d) for d in values(FLAT_CORR); init=0),
                sum(length(d) for d in values(FINE_PLANS); init=0),
                sum(length(d) for d in values(COARSE_PLANS); init=0),
                sum(length(d) for d in values(ANGULAR_SAME_PLANS); init=0),
                sum(length(d) for d in values(ANGULAR_REMOTE_PLANS); init=0)]
    st = Float64.(ANG_STATS)
    t = MPI.Allreduce(n, +, comm); s = MPI.Allreduce(st, +, comm); sm = MPI.Allreduce(st, MPI.MAX, comm)
    MPI.Comm_rank(comm) == 0 || return String[]
    return String["switches in effect: USE_FACE_PLANS=$(USE_FACE_PLANS[]) USE_LOCAL_NCF_PLANS=$(USE_LOCAL_NCF_PLANS[]) " *
        "USE_FLAT_FACES=$(USE_FLAT_FACES[]); plan entries built (sum over ranks): local fine $(Int(t[1])), " *
        "local coarse $(Int(t[2])), flat faces $(Int(t[3])), remote fine $(Int(t[4])), remote coarse $(Int(t[5]))",
        "angular plans: USE_ANGULAR_PLANS=$(USE_ANGULAR_PLANS[]); entries: same-rank $(Int(t[6])), cross-rank $(Int(t[7])); " *
        "measured phase: same-rank built $(Int(s[1])) evaluated $(Int(s[2])) (max rank $(Int(sm[2]))), " *
        "cross-rank built $(Int(s[3])) evaluated $(Int(s[4])) (max rank $(Int(sm[4]))), slot tables built $(Int(s[5]))"]
end

function bench_packed_cache(comm; npairs::Int=1000, n3::Int=125, reps_cached::Int=20000, reps_uncached::Int=3000)
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm); ni = node_info(comm)
    A0 = rand(n3, n3) + n3*I
    facs = [begin Ap = copy(A0); Ap[1,1] += 0.001*p; lu!(Ap) end for p in 1:npairs]     # 1000 x 125 KB per rank
    b = rand(n3); x = similar(b); W = similar(A0)
    for p in 1:20; x .= b; ldiv!(facs[p], x); end                                       # warm-up
    copyto!(W, A0); F = lu!(W); x .= b; ldiv!(F, x)
    MPI.Barrier(comm)
    t0 = time_ns()
    for r in 1:reps_cached
        p = 1 + (r*7919) % npairs                                                       # strided walk through the factors
        x .= b; ldiv!(facs[p], x)
    end
    tc = (time_ns() - t0) * 1e-9 / reps_cached
    MPI.Barrier(comm)
    t0 = time_ns()
    for r in 1:reps_uncached
        copyto!(W, A0); W[1,1] += 1e-3*r; F = lu!(W); x .= b; ldiv!(F, x)
    end
    tu = (time_ns() - t0) * 1e-9 / reps_uncached
    s = MPI.Allreduce(Float64[tc, tu], +, comm) ./ np
    m = MPI.Allreduce(Float64[tc, tu], MPI.MAX, comm)
    if rank == 0
        bytes = 8.0*n3*n3
        println("packed-node benchmark: $np ranks, $(ni.ranks_on_node) ranks on rank 0's node ($(ni.n_nodes) nodes): " *
            "cached ldiv! $(round(1e6*s[1]; digits=1)) us mean ($(round(1e6*m[1]; digits=1)) max); " *
            "uncached copy+lu!+ldiv! $(round(1e6*s[2]; digits=1)) us mean ($(round(1e6*m[2]; digits=1)) max); " *
            "uncached/cached $(round(s[2]/s[1]; digits=1))x; cached stream $(round(bytes/s[1]/1e9; digits=2)) GB/s per rank, " *
            "$(round(bytes/s[1]/1e9*ni.ranks_on_node; digits=1)) GB/s per node")
    end
    return nothing
end
