#=============================================================================
 schur_stage.jl -- the production stage solve, run through the SCALAR Schur
                   system instead of the five-field one.

 WHAT CHANGES, AND WHAT DOES NOT
 -------------------------------
 `imex3d_solve!` solves (I - gdt*A) U = B by GMRES over 5*Np unknowns. This
 solves the algebraically equivalent scalar system

     H[P] = rhs,      P = beta * Theta

 over Np unknowns, and rebuilds the five fields from P afterwards. Same answer:
 test/hevi/test_schur_solve.jl formed BOTH systems densely and solved them
 directly, and the five fields agreed to 3.3e-16 / 6.7e-14. Nothing about the
 time integrator, the splitting, the tolerance or the reference state changes;
 only the linear algebra the stage solve does.

 IT REQUIRES THE ADVECTIVE THETA ROW, AND THAT IS NOT A PREFERENCE
 -----------------------------------------------------------------
 With the FLUX-form Theta row, eliminating the momentum leaves a 2x2 system in
 (rho, P) -- 2*Np, not Np -- because Div[tb .* W m] carries a derivative of an
 unknown that continuity cannot cancel. The advective row is what makes
 `q := Theta - tb.*rho` pointwise and closes the reduction on one scalar. So
 `:imex_schur` forces `theta_advective` on the operator rather than offering it
 as a separate choice, and the deck cannot ask for one without the other.

 That changes f_imp, and therefore the SPLITTING -- but not the scheme's
 consistency or its answer to within the splitting error, because ark.jl never
 forms f_exp: it is `f(u) - f_imp(u)`, so whatever A is, the two halves still
 sum to the full right-hand side. test/hevi/test_theta_advective.jl measures the
 difference between the two Theta rows at 0.06% of the flux form.

 WHY THE COST FALLS
 ------------------
 Measured on the 64x64x60 deck (commit "Profile verdict"), 87.4% of a step is
 the stage solve, and of that: matvec 46.9%, banded solve 16.4%, gather/scatter
 2.2%, orthogonalisation 23.5%, MPI reduce 6.8%. Everything except the MPI
 reduce scales with the implicit field count, so 5 -> 1 removes most of it; the
 Krylov vectors go from 5*Np to Np, which is the part that also cuts memory
 traffic in the matvec. Step 5 measured the iteration count at the production
 anisotropy and stiffness at 3 against 8 for the five-field system, both
 preconditioned the way production can afford.
=============================================================================#

"""
    IMEX3DSchur

Everything the scalar stage solve needs, built once at setup.

`ws` is a SEPARATE GMRES workspace from the five-field one and is `npoin x 1`:
that is the whole point, and it is why this cannot share the existing
workspace. `DistributedInner` is field-count agnostic (it weights by node
multiplicity), so the SAME inner product object serves both.

`P` and `R` are `npoin`-vectors; `Pm` and `Rm` are `npoin x 1` RESHAPES of the
same memory, because `gmres_solve!` takes matrices while the Schur kernels take
vectors. A reshape shares storage, so there is no copy on either side of the
call.
"""
struct IMEX3DSchur{PCT}
    st::SchurState{Float64}
    pc::PCT                      # SchurColumnPrecond, or nothing
    ws::GMRESWorkspace
    P::Vector{Float64}
    R::Vector{Float64}
    Pm::Matrix{Float64}
    Rm::Matrix{Float64}
end

"""
    build_imex3d_schur(params, topo, comm, inputs, op, gdt, inner) -> IMEX3DSchur

Build the scalar workspace and, unless `:imex_precond` is `:none`, the one-field
column preconditioner.

The Krylov settings are read from the SAME deck keys the five-field solve uses
(`:imex_rtol`, `:imex_restart`, `:imex_maxiter`), so a deck that switches
`:imex_schur` on does not silently change its tolerance as well -- which would
make any before/after timing meaningless.
"""
function build_imex3d_schur(params, topo::ColumnTopology, comm::MPI.Comm, inputs,
                            op::HEVIOperator, gdt::Real, inner::DistributedInner;
                            lwall_flux::Bool = true)

    npoin = Int(params.mesh.npoin)

    op.theta_advective ||
        error("IMEX3D: the Schur stage solve needs the ADVECTIVE Theta row. ",
              "With the flux form the elimination leaves a 2x2 system in ",
              "(rho, P) rather than one scalar, so the reduction does not ",
              "close. build_imex3d sets this automatically when :imex_schur ",
              "is on; reaching this means the operator was built elsewhere.")

    pcmode  = Symbol(get(inputs, :imex_precond, :column))
    rtol    = Float64(get(inputs, :imex_rtol, 1.0e-8))
    restart = Int(get(inputs, :imex_restart, 20))
    maxiter = Int(get(inputs, :imex_maxiter, 200))

    # :imex_schur_kernel -- the bespoke scalar sweeps for H, DEFAULT ON.
    #
    # Off, `schur_H!` is two full five-field `hevi_apply_A!` calls. That is the
    # reference form and it is correct; it is also 2.05x one application, which
    # is why the first production profile of this path showed the matvec 46%
    # SLOWER than the five-field solve it replaced even while every other term
    # fell by 5-12x. On, the same H costs 0.36x one application -- 5.8x less --
    # and agrees with the reference to 1.9e-16 (test/hevi/test_schur_kernel.jl).
    #
    # The key exists so that a deck seeing a suspicious answer can go back to
    # the reference form WITHOUT editing source, and get a second, independent
    # statement of the same operator to compare against. That is the only thing
    # it is for; there is no accuracy reason to turn it off.
    usekern = Bool(get(inputs, :imex_schur_kernel, true))

    ws = GMRESWorkspace(npoin, 1, inner; m = restart, maxiter = maxiter,
                        rtol = rtol, atol = 1.0e-30)

    pc = pcmode === :column ?
         build_schur_column_precond(params, topo, comm, gdt;
                                    lwall_flux = lwall_flux) : nothing

    P = zeros(Float64, npoin); R = zeros(Float64, npoin)
    # The PRECONDITIONER deliberately gets no kernel: it runs on a vertical-only
    # operator, and the kernel is the full 3D form. `_schur_kernel` in schur.jl
    # turns that combination into an error rather than a quiet wrong band.
    kern = usekern ? SchurKernel(params) : nothing

    return IMEX3DSchur(SchurState(npoin, op.nimp; kern = kern), pc, ws, P, R,
                       reshape(P, npoin, 1), reshape(R, npoin, 1))
end

"""
    imex3d_solve_schur!(dst, src, params, gdt) -> dst

`dst = qe + (I - gdt*A)^-1 (src - qe)`, obtained through the scalar Schur
system rather than the five-field one.

THE ORDER OF THE THREE CALLS IS LOAD-BEARING. `schur_setup_rhs!` leaves `q_b`
and the known momentum `m_b` in the state, and `schur_recover!` reads both.
`schur_H!` in between writes only the transient `m`/`r`/`d` slots, never those
two -- which is what makes it safe to run the whole Krylov iteration between
setup and recovery. test/imex3d/test_schur_stage.jl asserts exactly that by
comparing against the five-field solve AFTER a full GMRES run, so a clobber
would show up as a wrong state rather than as a slow one.

EVERY RANK RUNS THE SAME NUMBER OF ITERATIONS, for the reason spelled out on
`imex3d_solve!`: the loop conditions come out of `MPI.Allreduce` and the
operator contains a halo exchange, so a rank-local branch here deadlocks.
"""
function imex3d_solve_schur!(dst, src, params, gdt::Real)

    imex  = params.imex
    sch   = imex.schur
    op    = imex.op
    npoin = Int(params.mesh.npoin)
    qe    = params.qp.qe
    g     = Float64(gdt)

    B = imex.B
    X = imex.X
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            B[ip, q] = src[off + ip] - qe[ip, ieq]
        end
    end

    # 5-field right-hand side -> scalar right-hand side, caching q_b and m_b.
    # ACCOUNTED SEPARATELY. This and the recovery below are per-SOLVE costs, not
    # per-iteration ones, so none of the matvec/precond/orthogonalise counters
    # can see them -- on the first production profile of this path they appeared
    # only as a 1.27 s/step hole between `sub-accounted` and the solve total,
    # 16% of the stage solve. They are also where an optimised scalar kernel for
    # H would show up second (setup_rhs is two full operator applications and
    # recover is one), which is the other reason to measure them apart.
    _tsr = hevi_tic()
    schur_setup_rhs!(sch.R, sch.st, B, params, op, g)
    if hevi_prof_on()
        HEVI_PROFILE.t_srhs += (time_ns() - _tsr) * 1e-9; HEVI_PROFILE.n_srhs += 1
    end

    st = sch.st
    matvec! = let params = params, op = op, g = g, st = st
        (W, V) -> begin
            schur_H!(vec(W), vec(V), params, op, g, st)
            return W
        end
    end
    precon! = sch.pc === nothing ? (Z -> Z) :
              let pc = sch.pc, params = params, g = g
                  Z -> (schur_precond!(vec(Z), pc, params, g); Z)
              end

    hevi_trace("    IMEX/Schur: entering GMRES on the scalar system, gdt=", g)
    iters, rel, converged = gmres_solve!(sch.Pm, sch.Rm, sch.ws, matvec!, precon!;
                                         warm = imex.warm_start)
    hevi_trace("    IMEX/Schur: GMRES done in ", iters, " iterations, rel. residual ", rel)

    if !converged && MPI.Comm_rank(imex.comm) == 0
        @warn string("IMEX3D/Schur: the scalar stage solve did not reach its ",
                     "tolerance -- ", iters, " iterations left a relative residual of ",
                     rel, " against :imex_rtol = ", sch.ws.rtol, ". Raise ",
                     ":imex_maxiter, raise :imex_restart, or lower Δt.") maxlog = 20
    end

    # rebuild the five fields from the pressure
    _trc = hevi_tic()
    schur_recover!(X, sch.P, sch.st, B, params, op, g)
    if hevi_prof_on()
        HEVI_PROFILE.t_srec += (time_ns() - _trc) * 1e-9; HEVI_PROFILE.n_srec += 1
    end

    dst === src || copyto!(dst, src)
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            dst[off + ip] = X[ip, q] + qe[ip, ieq]
        end
    end

    if IMEX_MONITOR[] && sch.ws.nsolve % IMEX_MONITOR_EVERY[] == 0 &&
       MPI.Comm_rank(imex.comm) == 0
        @printf(" [imex/schur] %d solves, %.1f Krylov iterations each on average; last %d at %.2e, worst %.2e\n",
                sch.ws.nsolve, sch.ws.niter / max(sch.ws.nsolve, 1), sch.ws.last_iters,
                sch.ws.last_relres, sch.ws.worst_relres)
        flush(stdout)
    end
    return dst
end
