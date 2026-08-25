#=============================================================================
 test_schur_precond.jl -- STEP 5 of the Schur complement work: the ONE-FIELD
 column preconditioner, and the like-for-like iteration count it unlocks.

 WHAT THIS IS FOR
 ----------------
 Step 4's sweep compared H against the full 5-field operator with "a column
 preconditioner" on both, and reported ~2x fewer iterations for H. That
 preconditioner was the exact inverse of the FULL 3D operator restricted to
 each column's nodes. Production does not have it: `imex3d_precond!` factorises
 the VERTICAL-ONLY operator. Block Jacobi by column keeps the diagonal part of
 the horizontal coupling; the vertical operator throws all of it away, so the
 two are different matrices and the step-4 ratio was measured against a
 preconditioner nobody can build at production cost.

 `build_schur_column_precond` builds the one that stands in the same relation
 to H as `imex3d_precond!` does to (I - lam*A). This file checks it is what it
 claims to be and then measures the ratio again, honestly.

 THE STIFFNESS IS SET FROM THE SPECTRUM, NOT PICKED
 --------------------------------------------------
 The first version of this file ran the measurement at the mock's default
 lam = 0.05 and reported "H 1 iteration, full 3". That is the SAME mistake
 step 2 recorded and fixed: at that lam the stage system is trivial
 (cond(I - lam*A) = 1.06) and the comparison measures nothing. Here lam is
 CHOSEN per mesh as `target/rate` with `rate = max|Im eig(A)|`, so every row of
 the sweep sits at a stated multiple of the acoustic stiffness and the
 production point (lam*rate = gamma*dt*c/h_z = 10.7 for LESICP2-64x64x60) is
 one of them.

 EVERYTHING DENSE IS SINGLE-RANK, AND THAT IS NOT A CONVENIENCE
 --------------------------------------------------------------
 `dense_scalar`/`dense_packed` loop to the RANK-LOCAL npoin, and every operator
 application inside them is collective. On 3 ranks the ranks therefore issue
 DIFFERENT numbers of Allreduces, and the job dies inside MPI with
 "Message truncated" pointing at setup_assembler -- nowhere near the cause.
 The dense assertions and the sweep are behind `if NR == 1` for that reason;
 the matrix-free round-trip at the end is what actually runs in parallel.

 ORDER MATTERS. The measurement is last on purpose: an iteration count is a
 plausible-looking number whatever the preconditioner actually contains, so the
 assertions that pin the matrix come first.

     julia --project=<env> test/hevi/test_schur_precond.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_schur_precond.jl
=============================================================================#

using Test, MPI, LinearAlgebra, Printf
MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "mock_sem.jl"))
for f in ("columns.jl", "vdiffusion.jl", "operator.jl", "factorize.jl",
          "schur.jl", "schur_precond.jl")
    include(joinpath(SRC, f))
end
say(a...) = RANK == 0 && (println(a...); flush(stdout))

# GMRES that counts iterations. Copied from test_schur_solve.jl together with
# the bug it had: the first version zeroed H[j+1,j] in the Givens rotation and
# then tested that same entry for breakdown, so it returned at j = 1 for every
# matrix including a cond-1e3 diagonal. `hnext` is kept separately for that
# reason -- see the commit that fixed it before believing any count.
function gmres_count(A, b; tol = 1e-4, maxit = 300)
    n = length(b); r = b - A*zeros(n); β = norm(r); β == 0 && return 0
    V = [r/β]; H = zeros(maxit+1, maxit); g = zeros(maxit+1); g[1] = β
    cs = zeros(maxit); sn = zeros(maxit)
    for j = 1:maxit
        w = A*V[j]
        for i = 1:j; H[i,j] = dot(V[i], w); w -= H[i,j]*V[i]; end
        hnext = norm(w); H[j+1,j] = hnext
        hnext > 1e-14 && push!(V, w/hnext)
        for i = 1:j-1
            t = cs[i]*H[i,j] + sn[i]*H[i+1,j]
            H[i+1,j] = -sn[i]*H[i,j] + cs[i]*H[i+1,j]; H[i,j] = t
        end
        ρ = hypot(H[j,j], H[j+1,j]); cs[j] = H[j,j]/ρ; sn[j] = H[j+1,j]/ρ
        H[j,j] = ρ; H[j+1,j] = 0.0
        g[j+1] = -sn[j]*g[j]; g[j] = cs[j]*g[j]
        abs(g[j+1])/β <= tol && return j
        hnext <= 1e-14 && return j      # PRE-rotation value; H[j+1,j] was zeroed
    end
    return maxit
end

"Dense matrix of a scalar operator `f!(out, in)`, column by column."
function dense_scalar(f!, n)
    M = zeros(n, n); e = zeros(n)
    for j = 1:n
        fill!(e, 0.0); e[j] = 1.0
        M[:, j] = f!(zeros(n), e)
    end
    return M
end

"Dense matrix of a packed npoin x nimp operator `f!(W, V)`."
function dense_packed(f!, n, nimp)
    m = n * nimp; M = zeros(m, m)
    V = zeros(n, nimp); W = zeros(n, nimp)
    for j = 1:m
        fill!(V, 0.0); V[j] = 1.0; fill!(W, 0.0); f!(W, V)
        @inbounds for i = 1:m; M[i, j] = W[i]; end
    end
    return M
end

#-----------------------------------------------------------------------------
# One mesh: build everything, and return the pieces both the assertions and the
# sweep need. `lam` is set later, from the spectrum.
#-----------------------------------------------------------------------------
struct Case
    params; topo; op; st; N::Int; ngl::Int
    Araw::Matrix{Float64}       # the raw operator A, 5 fields, dense
    rate::Float64               # max |Im eig(A)| -- the acoustic stiffness
end

function build_case(; nelx, nely, nelz, Lx, Ly, Lz, p)
    params = build_mock_params(nelx=nelx, nely=nely, nelz=nelz, p=p,
                               Lx=Lx, Ly=Ly, Lz=Lz, comm=COMM, base=:stratified)
    topo = build_column_topology(params.mesh, COMM)
    # ADVECTIVE: the reduction closes on one scalar only with this row.
    op = build_hevi_operator(params, topo, [1,2,3,4,5]; lwall_flux=true,
                             full=true, theta_advective=true)
    N  = Int(params.mesh.npoin)
    st = SchurState(N, op.nimp)
    # DENSE ONLY ON ONE RANK. `dense_packed` loops to the RANK-LOCAL npoin and
    # every hevi_apply_A! inside it is collective, so on >1 rank the ranks issue
    # different numbers of Allreduces and the job dies inside MPI with a
    # truncated-message error nowhere near the cause. Everything that consumes
    # Araw is already behind `if NR == 1`.
    NR > 1 && return Case(params, topo, op, st, N, Int(params.mesh.ngl),
                          zeros(0, 0), NaN)
    Araw = dense_packed((W, V) -> hevi_apply_A!(W, V, params, op), N, op.nimp)
    rate = maximum(abs, imag.(eigvals(Araw)))
    return Case(params, topo, op, st, N, Int(params.mesh.ngl), Araw, rate)
end

"""
Iteration counts for both operators at stiffness `lam`, each preconditioned by
the inverse of its OWN vertical operator factorised per column -- which is what
production can afford for either one.
"""
function counts_at(c::Case, lam::Float64)
    params, topo, op, st, N = c.params, c.topo, c.op, c.st, c.N

    # --- the scalar system: H, preconditioned by the vertical H
    pc = build_schur_column_precond(params, topo, COMM, lam)
    Hfull  = dense_scalar((o, e) -> schur_H!(o, e, params, op, lam, st), N)
    Mschur = dense_scalar((o, e) -> (copyto!(o, e);
                                     schur_precond!(o, pc, params, lam)), N)

    # --- the 5-field system: (I - lam*A), preconditioned by the vertical A,
    #     applied exactly as imex3d_precond! applies it.
    pvars = [1,2,3,4,5]
    opv5  = build_hevi_operator(params, topo, pvars; lwall_flux=true,
                               full=false, theta_advective=true)
    owner, own = assign_column_owners(topo, COMM)
    cc5  = build_column_comm(topo, owner, own, COMM, length(pvars))
    fac5 = build_column_factorization(params, opv5, cc5, topo, lam)
    m = N * op.nimp
    Afull = Matrix(1.0I, m, m) .- lam .* c.Araw
    M5 = dense_packed((W, V) -> begin
             copyto!(W, V)
             column_gather!(cc5, W, pvars)
             for ic = 1:cc5.nown
                 LAPACK.gbtrs!('N', fac5.kl, fac5.ku, fac5.n, fac5.AB[ic],
                               fac5.ipiv[ic], view(cc5.X, :, ic))
             end
             column_scatter!(cc5, W, pvars)
         end, N, op.nimp)

    gip = c.params.mesh.ip2gip
    bH = [sinpi(1e-2*gip[i]) for i = 1:N]
    bF = vec([sinpi(1e-2*(gip[i] + 7j)) for i = 1:N, j = 1:op.nimp])
    return gmres_count(Hfull * Mschur, bH), gmres_count(Afull * M5, bF)
end

#-----------------------------------------------------------------------------
const BASE = build_case(nelx=2, nely=2, nelz=5, Lx=400.0, Ly=400.0, Lz=1000.0, p=2)
# Production stiffness: gamma*dt*c/h_z = 10.7 on LESICP2-64x64x60. Written as a
# LITERAL rather than as 10.7/BASE.rate because `rate` needs a dense
# eigendecomposition, which only one rank can compute -- and a lam that differed
# with the rank count would make the parallel check below meaningless. The
# 1-rank testset asserts the literal still matches the spectrum, so it cannot
# drift if the mock changes.
const LAM = 1.5039
const LAM_TARGET = 10.7
const N   = BASE.N
const NGL = BASE.ngl

# The band, captured UNFACTORISED through the hook -- the only moment it can be
# compared against the operator, since gbtrf! overwrites it.
band = Ref{Any}(nothing)
pc = build_schur_column_precond(BASE.params, BASE.topo, COMM, LAM;
                                verify_hook = (AB, kl, ku, n) ->
                                    (band[] = (deepcopy(AB), kl, ku, n)))

say("\n=== Schur step 5: the one-field column preconditioner, $NR rank(s) ===")
NR == 1 ?
    say(@sprintf("  mock stiffness rate = %.4g, so lam = %.4g puts lam*rate at %.4g",
                 BASE.rate, LAM, LAM*BASE.rate)) :
    say(@sprintf("  lam = %.4g; the spectrum needs a dense eigensolve, so it is",
                 LAM), " measured and asserted on 1 rank only")

@testset "Schur column preconditioner ($NR ranks)" begin

    params, topo, op, st = BASE.params, BASE.topo, BASE.op, BASE.st

if NR == 1
    # The vertical-only Schur operator, densely. This is what the
    # preconditioner is supposed to invert. INSIDE the guard: dense_scalar
    # loops to the rank-local npoin and schur_H! is collective, so building
    # this on every rank deadlocks the moment the ranks own different numbers
    # of nodes -- which is exactly how the first 3-rank run of this file hung.
    Hv = dense_scalar((o, e) -> schur_H!(o, e, params, pc.opv, LAM, pc.st), N)

    lev = level_of_node(topo, N)
    colof = zeros(Int, N)
    for ic = 1:topo.ncol, il = 1:topo.nlev
        ip = topo.node[il, ic]; ip != 0 && (colof[ip] = ic)
    end

    @testset "lam really is at production stiffness" begin
        # LAM is a literal so it cannot vary with the rank count. This is what
        # stops it drifting: if the mock's spectrum ever changes, the literal
        # stops meaning "production stiffness" and this fails rather than
        # quietly moving the whole sweep to a different regime -- which is the
        # failure that made the first version of this file report
        # "H 1 iteration, full 3" and measure nothing.
        say(@sprintf("  rate = %.6g, lam*rate = %.4g (target %.4g)",
                     BASE.rate, LAM * BASE.rate, LAM_TARGET))
        @test isapprox(LAM * BASE.rate, LAM_TARGET; rtol = 1e-3)
    end

    @testset "the vertical Schur operator is block-diagonal by column" begin
        # The PREMISE of a column preconditioner. If the vertical operator
        # coupled columns, factorising it per column would be an approximation
        # of an approximation and the band would not be the operator at all.
        off = 0.0
        for i = 1:N, j = 1:N
            colof[i] == colof[j] || (off = max(off, abs(Hv[i,j])))
        end
        scale = maximum(abs, Hv)
        say(@sprintf("  off-column coupling in H_v: %.3e (relative %.3e)",
                     off, off/scale))
        @test off < 1e-12 * scale
    end

    @testset "the stencil reaches exactly 2(ngl-1) levels" begin
        # THE ONE THING THAT ALIASES SILENTLY. The colouring period is
        # 4(ngl-1)+1 because H = Div . A^-1 . Grad reaches ngl-1 twice. If the
        # true reach were larger, the probe would fold outer entries onto the
        # wrong columns and still produce a working, wrong preconditioner.
        w = schur_column_reach(NGL)
        @test w == 2 * (NGL - 1)
        beyond = 0.0; at_edge = 0.0
        for i = 1:N, j = 1:N
            colof[i] == colof[j] || continue
            d = abs(lev[i] - lev[j])
            d >  w && (beyond  = max(beyond,  abs(Hv[i,j])))
            d == w && (at_edge = max(at_edge, abs(Hv[i,j])))
        end
        scale = maximum(abs, Hv)
        say(@sprintf("  |H_v| beyond %d levels: %.3e; at exactly %d: %.3e (scale %.3e)",
                     w, beyond, w, at_edge, scale))
        @test beyond < 1e-12 * scale
        # ... and the outermost diagonal is NOT empty, or the claim would be
        # vacuous and a narrower band would have done.
        @test at_edge > 1e-9 * scale
    end

    @testset "the assembled band IS the vertical operator" begin
        # `assemble_schur_column_band!` recovers H_v by probing. Compare every
        # stored entry against the dense operator: this is the check that the
        # colouring, the band indexing and the zero seed are all right at once.
        AB, kl, ku, n = band[]
        diagrow = kl + ku + 1
        worst = 0.0
        for ic = 1:pc.cc.nown
            gc = pc.cc.own_gcol[ic]
            for il = 1:n, jl = max(1, il-kl):min(n, il+ku)
                ip = topo.node[il, gc]; jp = topo.node[jl, gc]
                (ip == 0 || jp == 0) && continue
                worst = max(worst, abs(AB[ic][diagrow + il - jl, jl] - Hv[ip, jp]))
            end
        end
        say(@sprintf("  max |assembled band - dense H_v| = %.3e", worst))
        @test worst < 1e-10 * maximum(abs, Hv)
    end

    @testset "the preconditioner actually inverts the vertical operator" begin
        # ||I - H_v * M|| must be ~0: M is supposed to be H_v^-1 exactly, since
        # H_v is block-diagonal by column and each block is factorised whole.
        M = dense_scalar((o, e) -> (copyto!(o, e);
                                    schur_precond!(o, pc, params, LAM)), N)
        err = opnorm(Matrix(1.0I, N, N) .- Hv * M, Inf)
        say(@sprintf("  ||I - H_v * M||_inf = %.3e", err))
        @test err < 1e-8
    end

    @testset "the exact inverse converges in one iteration" begin
        # Validates the PRECONDITIONED PATH itself before any count taken
        # through it is believed: feeding GMRES the exact inverse must give 1.
        Hfull = dense_scalar((o, e) -> schur_H!(o, e, params, op, LAM, st), N)
        gip = params.mesh.ip2gip
        @test gmres_count(Hfull * inv(Hfull), [sinpi(1e-2*gip[i]) for i = 1:N]) == 1
    end

    @testset "why this differs from step 4: the two preconditioners compared" begin
        # Step 4 measured the 1:1 mesh at 47 iterations for H and called it a
        # property of the ISOTROPIC MESH ("a column preconditioner is the wrong
        # tool there"). It is not. It is a property of the PRECONDITIONER step 4
        # used -- the exact inverse of the FULL operator restricted to a
        # column's nodes, i.e. block Jacobi, which production cannot build
        # anyway. Both are applied to the SAME operator here so the comparison
        # cannot be confounded by setup differences.
        c = build_case(nelx=2, nely=2, nelz=5, Lx=400.0, Ly=400.0, Lz=1000.0, p=2)
        lam = 2.0 / c.rate
        pcv = build_schur_column_precond(c.params, c.topo, COMM, lam)
        Hf  = dense_scalar((o, e) -> schur_H!(o, e, c.params, c.op, lam, c.st), c.N)
        Mv  = dense_scalar((o, e) -> (copyto!(o, e);
                                      schur_precond!(o, pcv, c.params, lam)), c.N)
        CB  = let M = zeros(size(Hf))
            for ic = 1:c.topo.ncol
                idx = [c.topo.node[il, ic] for il = 1:c.topo.nlev
                       if c.topo.node[il, ic] != 0]
                isempty(idx) || (M[idx, idx] = inv(Hf[idx, idx]))
            end
            M
        end
        Id = Matrix(1.0I, c.N, c.N)
        ev = opnorm(Id .- Hf*Mv, Inf)
        eb = opnorm(Id .- Hf*CB, Inf)
        gip = c.params.mesh.ip2gip
        bH  = [sinpi(1e-2*gip[i]) for i = 1:c.N]
        iv, ib = gmres_count(Hf*Mv, bH), gmres_count(Hf*CB, bH)
        say(@sprintf("  1:1 mesh, lam*rate = 2: vertical-operator precond  ||I-HM|| %.3g, %d its",
                     ev, iv))
        say(@sprintf("                          block-Jacobi-of-full       ||I-HM|| %.3g, %d its",
                     eb, ib))
        # The vertical operator is a consistent sub-operator of H; block Jacobi
        # truncates the coupling instead, and on an isotropic mesh what it
        # truncates is as large as what it keeps.
        @test iv <= ib
    end

    @testset "THE MEASUREMENT: H vs full, both preconditioned as production is" begin
        # Swept over BOTH parameters that govern this: the acoustic stiffness
        # lam*rate, and the grid anisotropy h_x/h_z that decides how much a
        # COLUMN preconditioner can possibly do for either operator.
        #
        # Rows are (elements in m) x (lam*rate). The production point is the
        # 4:1 mesh at lam*rate = 10.7.
        meshes = (("  1:1  200x200x200", (2, 2, 5, 400.0, 400.0, 1000.0)),
                  ("  4:1  400x400x100", (2, 2, 10, 800.0, 800.0, 1000.0)),
                  (" 16:1  800x800x 50", (2, 2, 20, 1600.0, 1600.0, 1000.0)))
        say("\n  anisotropy          lam*rate      H   full   ratio")
        worst_ratio = Inf
        for (name, (nx, ny, nz, Lx, Ly, Lz)) in meshes
            c = build_case(nelx=nx, nely=ny, nelz=nz, Lx=Lx, Ly=Ly, Lz=Lz, p=2)
            for target in (2.0, 10.7)
                lam = target / c.rate
                iH, iF = counts_at(c, lam)
                say(@sprintf("%s %8.2f   %4d   %4d   %5.2fx",
                             name, target, iH, iF, iF/max(iH,1)))
                worst_ratio = min(worst_ratio, iF/max(iH,1))
                # The claim under test is only that the scalar system is no
                # HARDER than the 5-field one once both are preconditioned the
                # way production can afford. The SIZE of the win is reported,
                # not asserted -- the mock is a few hundred dof at p=2 against
                # 15.9M at p=4, and the iteration ratio is the part least
                # likely to transfer.
                @test iH <= iF
            end
        end
        say(@sprintf("\n  worst ratio over the sweep: %.2fx", worst_ratio))
    end
end # NR == 1

    @testset "the preconditioner inverts H_v at ANY rank count" begin
        # Matrix-free, so it runs identically on 1 and on 3 ranks -- the dense
        # checks above cannot, because they loop to the RANK-LOCAL npoin while
        # every operator application inside is collective.
        #
        # H_v is block-diagonal by column and each block is factorised whole,
        # so M is its EXACT inverse: apply H_v to a field, solve, get the field
        # back. A split column is assembled from pieces on different ranks and
        # solved on its owner, so this failing on 3 ranks and passing on 1 is
        # precisely the gather/scatter bug it exists to catch.
        gip = params.mesh.ip2gip
        P0 = [sinpi(1e-2 * gip[i]) + 0.3cospi(3e-2 * gip[i]) for i = 1:N]
        Q  = schur_H!(zeros(N), P0, params, pc.opv, LAM, pc.st)
        R  = schur_precond!(copy(Q), pc, params, LAM)
        err = maximum(abs, R .- P0)
        tot = MPI.Allreduce(err, MPI.MAX, COMM)
        say(@sprintf("  max |M[H_v[P]] - P| over all ranks: %.3e", tot))
        @test tot < 1e-9 * maximum(abs, P0)
    end
end

MPI.Barrier(COMM)
say("=== step 5 done on $NR rank(s) ===\n")
