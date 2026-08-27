#=  End-to-end: the PRODUCTION scalar stage solve with the kernel on and off.

    This is not a projection. The kernel and the reference form are the same
    operator to 1.9e-16, so both arms run the SAME number of Krylov iterations
    on the same right-hand side to the same tolerance -- the only difference is
    how each matvec, each setup_rhs and each recover is computed. The ratio
    therefore transfers directly to any mesh and any rank count, unlike an
    iteration-count comparison, which does not (the cold-start self-check on
    this path inverted the warm-start production truth once already).

    Measured, 1 rank, p = 4, npoin = 14161:

      iterations   8 and 8            (they must match, and are asserted to)
      states       2.952e-20 apart
      reference H  0.10806 s per stage solve
      kernel H     0.02664 s per stage solve
      speedup      4.06 x on the WHOLE stage solve

    4.06x rather than the 5.8x the H matvec alone gains, because the stage
    solve also contains the preconditioner, the orthogonalisation and the MPI
    reduce, none of which the kernel touches. The step will gain less again:
    rhs!, the explicit stages and the DSS outside the solve are unchanged.

      julia --project=. tools/schur_kernel_e2e.jl
=#
using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq
MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const SRC = joinpath(@__DIR__, "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "..", "test", "hevi", "mock_sem.jl"))
for f in ("columns.jl", "vdiffusion.jl", "operator.jl", "factorize.jl", "acoustic.jl",
          "ark.jl", "hevi.jl", "krylov.jl", "precond_api.jl",
          "schur.jl", "schur_precond.jl", "schur_stage.jl", "imex3d.jl")
    include(joinpath(SRC, f))
end
say(a...) = RANK == 0 && (println(a...); flush(stdout))

const NELX, NELY, NELZ, P = 4, 4, 12, 4
const GDT = 0.30
params0 = build_mock_params(nelx=NELX, nely=NELY, nelz=NELZ, p=P,
                            Lx=4000.0, Ly=4000.0, Lz=6000.0, comm=COMM,
                            base=:stratified)
topo  = build_column_topology(params0.mesh, COMM)
npoin = Int(params0.mesh.npoin)
op = build_hevi_fast_operator(params0, topo; lwall_flux = true, theta_advective = true)
inner = build_distributed_inner(op.dss_cache, npoin, COMM)

const RTOL = 1.0e-8
mkin(k) = Dict{Symbol,Any}(:imex_rtol => RTOL, :imex_restart => 20,
                           :imex_maxiter => 200, :imex_precond => :column,
                           :imex_schur_kernel => k)

sch_k = build_imex3d_schur(params0, topo, COMM, mkin(true),  op, GDT, inner)
sch_r = build_imex3d_schur(params0, topo, COMM, mkin(false), op, GDT, inner)
@assert sch_k.st.w.kern !== nothing && sch_r.st.w.kern === nothing

ws = GMRESWorkspace(npoin, op.nimp, inner; m = 20, maxiter = 200, rtol = RTOL, atol = 1e-30)
mk(s) = IMEX3DCache(topo, op, nothing, s, ws,
                    zeros(Float64, npoin, op.nimp), zeros(Float64, npoin, op.nimp),
                    imex3d_fimp!, imex3d_solve_schur!, :ARS343, GDT, :RS, 1, false, COMM)
pk = merge(params0, (imex = mk(sch_k),))
pr = merge(params0, (imex = mk(sch_r),))

gip = params0.mesh.ip2gip
src = zeros(Float64, npoin * 5)
for (q, ieq) in enumerate(op.vars), ip = 1:npoin
    src[(ieq-1)*npoin + ip] = params0.qp.qe[ip, ieq] +
        1.0e-3 * sinpi(3.0e-3 * (gip[ip] + 7q)) * (q == 1 ? 1.0 : 0.3)
end
dk = similar(src); dr = similar(src)

function bench(f, nb = 5, n = 8)
    f(); m = Inf
    for _ = 1:nb
        t = time_ns(); for _ = 1:n; f(); end
        m = min(m, (time_ns() - t) * 1e-9 / n)
    end
    m
end

imex3d_solve_schur!(dk, src, pk, GDT); ik = pk.imex.schur.ws.last_iters
imex3d_solve_schur!(dr, src, pr, GDT); ir = pr.imex.schur.ws.last_iters
say(@sprintf("npoin %d  iterations: kernel %d, reference %d  (must match)", npoin, ik, ir))
say(@sprintf("states agree to %.3e (relative)",
             maximum(abs, dk .- dr) / max(maximum(abs, dr), eps())))

tk = bench(() -> imex3d_solve_schur!(dk, src, pk, GDT))
tr = bench(() -> imex3d_solve_schur!(dr, src, pr, GDT))
say(@sprintf("\nFULL STAGE SOLVE, %d iterations each:", ik))
say(@sprintf("  reference H : %.5f s", tr))
say(@sprintf("  kernel H    : %.5f s", tk))
say(@sprintf("  speedup     : %.2f x  on the whole stage solve", tr/tk))
