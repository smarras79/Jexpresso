# RT timing record: per-call phase times, solve statistics, problem size and memory for get_rt_heat_fluxes_sweeps!.
#
# Purpose: compare the cost of RT configurations (uniform vs adaptive spatial/angular meshes, layouts) from measured
# numbers rather than job wall time. On by default (:RT_timing => true); the cost is a few Allreduces per RT call.
#
# Every RT call prints one block on rank 0 and appends rows to <output_dir>/profile/rt_timing.csv in long format:
#   call, nranks, quantity, max, mean, min        (max/mean/min over ranks; for global sums all three are the sum)
# Quantities:
#   phase_<name>        wall time of each phase of the call [s] (infra, atmos, sw_optics, adapt, sw_case, sw_solve, ...)
#   total               wall time of the whole call [s]
#   nelem, ang_dofs     global spatial elements, global Σ_elements (angular nodes); unknowns = ang_dofs * ngl^3
#   unknowns_rank       unknowns per rank (max/mean/min: partition balance)
#   <band>_matvecs, <band>_setup_s, <band>_solve_s, <band>_matvec1_s, <band>_steady_matvec_s,
#   <band>_steady_src_s, <band>_steady_sweep_s   (steady = matvecs 3..n, mean; per-matvec times are max over ranks)
#   <band>_sweep_rank_s total steady sweep time per rank (max/min: sweep load balance)
#   <band>_cache_mb, <band>_rss_gb               element LU cache and peak RSS (max over ranks)
#   rss_gb, live_gb     peak RSS and live Julia heap at the end of the call [GB]

using Printf

const RT_TIMING_NCALLS = Ref(0)
# All rows written so far (rank 0). The file is rewritten in full on every call: the first RT call happens during
# the integrator warm-up, and the output directory can be cleared after it, which would otherwise lose the cold solve.
const RT_TIMING_ROWS = String[]

mutable struct RTTiming
    names::Vector{Symbol}
    secs::Vector{Float64}
    solves::Vector{Pair{Symbol,Dict{Symbol,Any}}}
    t_start::UInt64
    t_lap::UInt64
end
RTTiming() = (t = time_ns(); RTTiming(Symbol[], Float64[], Pair{Symbol,Dict{Symbol,Any}}[], t, t))

# Close the current phase under `name` and start the next one.
function rt_lap!(tm::RTTiming, name::Symbol)
    now = time_ns()
    push!(tm.names, name); push!(tm.secs, (now - tm.t_lap) * 1e-9)
    tm.t_lap = now
    return nothing
end

# A Dict to pass as run_unified_case_solve(...; report=...) for band `band` (:sw, :lw).
function rt_solve_report!(tm::RTTiming, band::Symbol)
    d = Dict{Symbol,Any}()
    push!(tm.solves, band => d)
    return d
end

rt_timing_enabled(inputs) = get(inputs, :RT_timing, true) == true

# Collective: every rank must call it with the same phase names and solves.
function rt_timing_report!(tm::RTTiming, mesh, infra, comm, inputs)
    RT_TIMING_NCALLS[] += 1
    call = RT_TIMING_NCALLS[]
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    total = (time_ns() - tm.t_start) * 1e-9

    rows = Tuple{String,Float64,Float64,Float64}[]          # (quantity, max, mean, min)

    # phases + total: max / mean / min over ranks
    v = vcat(tm.secs, total)
    vmax = MPI.Allreduce(v, MPI.MAX, comm); vsum = MPI.Allreduce(v, MPI.SUM, comm); vmin = MPI.Allreduce(v, MPI.MIN, comm)
    for (i, n) in enumerate(tm.names)
        push!(rows, ("phase_$(n)", vmax[i], vsum[i]/np, vmin[i]))
    end
    push!(rows, ("total", vmax[end], vsum[end]/np, vmin[end]))

    # problem size
    ngl = mesh.ngl
    ndir_loc = mesh.nelem == 0 ? 0 : sum(infra.n_dirs[1:mesh.nelem])
    unk_loc = Float64(ndir_loc) * ngl^3
    g = MPI.Allreduce(Float64[mesh.nelem, ndir_loc, unk_loc], MPI.SUM, comm)
    umax = MPI.Allreduce(unk_loc, MPI.MAX, comm); umin = MPI.Allreduce(unk_loc, MPI.MIN, comm)
    push!(rows, ("nelem", g[1], g[1], g[1]))
    push!(rows, ("ang_dofs", g[2], g[2], g[2]))
    push!(rows, ("unknowns", g[3], g[3], g[3]))
    push!(rows, ("unknowns_rank", umax, g[3]/np, umin))

    # solves (already reduced inside run_unified_case_solve)
    for (band, d) in tm.solves
        haskey(d, :nmv) || continue
        nmv = d[:nmv]; st = nmv >= 3 ? (3:nmv) : (1:nmv)
        mt = d[:matvec_times]; sr = d[:src_times]; sw = d[:sweep_times]
        b = String(band)
        put(q, x) = push!(rows, ("$(b)_$(q)", x, x, x))
        put("matvecs", nmv)
        put("setup_s", MPI.Allreduce(d[:setup_time], MPI.MAX, comm))
        put("solve_s", d[:solve_time])
        put("matvec1_s", nmv >= 1 ? mt[1] : 0.0)
        put("steady_matvec_s", nmv >= 1 ? sum(mt[st])/length(st) : 0.0)
        put("steady_src_s", nmv >= 1 ? sum(sr[st])/length(st) : 0.0)
        put("steady_sweep_s", nmv >= 1 ? sum(sw[st])/length(st) : 0.0)
        push!(rows, ("$(b)_sweep_rank_s", d[:sweep_rank_max], NaN, d[:sweep_rank_min]))
        put("cache_mb", d[:cache_mb])
        put("rss_gb", d[:rss_max_gb])
    end

    # memory at the end of the call
    m = Float64[Sys.maxrss()/1e9, Base.gc_live_bytes()/1e9]
    mmax = MPI.Allreduce(m, MPI.MAX, comm); msum = MPI.Allreduce(m, MPI.SUM, comm); mmin = MPI.Allreduce(m, MPI.MIN, comm)
    push!(rows, ("rss_gb", mmax[1], msum[1]/np, mmin[1]))
    push!(rows, ("live_gb", mmax[2], msum[2]/np, mmin[2]))

    rank == 0 || return rows

    # log block
    f(x) = isnan(x) ? "-" : (abs(x) >= 1e5 ? @sprintf("%.4g", x) : @sprintf("%.3f", x))
    println("RT timing (call $call, $np ranks): total $(f(vmax[end])) s; " *
            "$(Int(g[1])) elements, $(Int(g[2])) angular dofs, $(Int(g[3])) unknowns " *
            "(per rank max/min $(Int(umax))/$(Int(umin)))")
    println("  phase [s]      max      mean       min")
    for (i, n) in enumerate(tm.names)
        @printf("  %-10s %9.3f %9.3f %9.3f\n", n, vmax[i], vsum[i]/np, vmin[i])
    end
    for (band, d) in tm.solves
        haskey(d, :nmv) || continue
        nmv = d[:nmv]; st = nmv >= 3 ? (3:nmv) : (1:nmv); mt = d[:matvec_times]
        imb = d[:sweep_rank_min] > 0 ? d[:sweep_rank_max] / d[:sweep_rank_min] : NaN
        println("  $(band): $nmv matvecs, setup $(f(rows_lookup(rows, "$(band)_setup_s"))) s, " *
                "solve $(f(d[:solve_time])) s, matvec 1 $(f(nmv >= 1 ? mt[1] : 0.0)) s, " *
                "steady matvec $(f(nmv >= 1 ? sum(mt[st])/length(st) : 0.0)) s, sweep imbalance max/min $(f(imb)), " *
                "cache $(f(d[:cache_mb])) MB, peak RSS $(f(d[:rss_max_gb])) GB")
    end
    println("  memory at end: peak RSS max $(f(mmax[1])) GB, live heap max $(f(mmax[2])) GB")
    flush(stdout)

    # CSV (long format)
    od = get(inputs, :output_dir, "none")
    if od isa AbstractString && od != "none" && !isempty(od)
        try
            dir = joinpath(od, "profile"); mkpath(dir)
            path = joinpath(dir, "rt_timing.csv")
            for (q, a, b, c) in rows
                push!(RT_TIMING_ROWS, @sprintf("%d,%d,%s,%.6g,%.6g,%.6g", call, np, q, a, b, c))
            end
            open(path, "w") do io
                println(io, "call,nranks,quantity,max,mean,min")
                foreach(r -> println(io, r), RT_TIMING_ROWS)
            end
        catch err
            @warn "RT timing: could not write rt_timing.csv" exception=err
        end
    end
    return rows
end

rows_lookup(rows, q) = (i = findfirst(r -> r[1] == q, rows); i === nothing ? NaN : rows[i][2])
