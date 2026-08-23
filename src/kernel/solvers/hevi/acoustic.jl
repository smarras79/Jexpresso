#=============================================================================
 acoustic.jl -- the FAST operator of an acoustic-substepping split, and its
 verification.

 WHAT THIS IS FOR
 ----------------
 HEVI removes the vertical acoustic term from the explicit eigenvalue budget
 and nothing else, so its gain is bounded by the grid's acoustic anisotropy.
 Getting past that means removing the HORIZONTAL acoustic term too. A 3D
 implicit solve would do it and is not affordable; acoustic substepping does
 it by integrating the acoustic subsystem with its own small step inside an
 outer step that then sees no sound at all.

 The piece both approaches need is an operator that IS the acoustic subsystem:

     f_fast(u) = A_full (u - qe)

 where A_full is the complete linear acoustic-gravity operator in all three
 directions. It is the same object as the HEVI operator with the ξ and η
 derivative sweeps put back -- same linearised fluxes, same reference state,
 same floor/lid treatment, same DSS and inverse mass matrix -- so it is built
 by `build_hevi_operator(...; full = true)` rather than written twice.

 The slow part is then, exactly as in HEVI,

     f_slow(u) = rhs!(u) - f_fast(u)

 computed by subtraction, so no source term, boundary condition, SGS closure,
 wall model or sponge can be lost or double-counted in the split.

 WHY IT IS VERIFIED THE WAY IT IS
 --------------------------------
 A fast operator that is not the exact acoustic part of the explicit
 discretisation does not fail loudly. The split stays consistent whatever the
 fast operator is (f_slow absorbs the difference), so a convergence test still
 passes; what happens instead is that a residual fast mode stays in the slow
 part and the outer step stays pinned by sound, with nothing in the output to
 say why. The check below is therefore against the operator the code already
 trusts, not against a paper derivation.
=============================================================================#

"""
    build_hevi_fast_operator(params, topo; lwall_flux = true) -> HEVIOperator

The complete linear acoustic-gravity operator, for use as the fast operator of
a substepping split. Carries all five equations, because the horizontal
pressure gradient and horizontal mass flux live in the ones HEVI leaves out.

Never factorised: it is not column-local, and in a substepping scheme it is
only ever applied.
"""
function build_hevi_fast_operator(params, topo::ColumnTopology; lwall_flux::Bool = true)
    return build_hevi_operator(params, topo, [1, 2, 3, 4, 5];
                               lwall_flux = lwall_flux, full = true)
end

"""
    hevi_fast_rhs!(du, u, params, op) -> du

`f_fast` in the flat `du` layout the rest of the code uses: the full acoustic
operator applied to the deviation `u - qe`.

Identical in structure to `hevi_implicit_rhs!`; kept separate so the two can
be used in the same run (a substepping scheme wants the fast operator for the
substeps and the vertical one for the implicit half of each substep).
"""
function hevi_fast_rhs!(du, u, params, op::HEVIOperator)

    npoin = Int(params.mesh.npoin)
    qe    = params.qp.qe

    V = op.V
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            V[ip, q] = u[off + ip] - qe[ip, ieq]
        end
    end

    hevi_apply_A!(op.R, V, params, op)

    fill!(du, 0.0)
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            du[off + ip] = op.R[ip, q]
        end
    end
    return du
end

"""
    hevi_verify_fast(params, opfull, opvert, topo) -> NamedTuple

Check the full acoustic operator against the vertical one the code already
trusts.

THE TEST. Give both operators a deviation field that depends on the vertical
level ONLY -- exactly uniform in x and y, so its horizontal derivatives vanish
discretely, not just to truncation error. On a mesh whose ζ lines are vertical
the ξ and η sweeps then contribute nothing, and the full operator must
reproduce the vertical operator to round-off in (ρ, ρw, ρθ).

That is a real test of the new sweeps: getting a metric term, a flux slot or a
derivative index wrong in them shows up as a mismatch here, because those
terms are supposed to cancel to zero and a wrong one does not.

It also reports what the horizontal terms do on a field that varies in x,
which must be NON-zero in (ρu, ρv) -- a full operator that silently reduced to
the vertical one would pass the first test and fail this one.

Returns `(rel_vertical, horiz_response, ok)`.
"""
function hevi_verify_fast(params, opfull::HEVIOperator, opvert::HEVIOperator,
                          topo::ColumnTopology)

    npoin = Int(params.mesh.npoin)
    nlev  = topo.nlev

    # ---- 1. horizontally uniform field: full must equal vertical ------------
    Vf = zeros(Float64, npoin, opfull.nimp)
    Vv = zeros(Float64, npoin, opvert.nimp)
    @inbounds for ic = 1:topo.ncol, il = 1:nlev
        ip = topo.node[il, ic]
        ip == 0 && continue
        a = sinpi(2.0 * (il - 1) / max(nlev - 1, 1))
        b = cospi(2.0 * (il - 1) / max(nlev - 1, 1))
        for (q, ieq) in enumerate(opfull.vars)
            Vf[ip, q] = ieq == 1 ? a : (ieq == 4 ? b : (ieq == 5 ? a * b : 0.0))
        end
        for (q, ieq) in enumerate(opvert.vars)
            Vv[ip, q] = ieq == 1 ? a : (ieq == 4 ? b : (ieq == 5 ? a * b : 0.0))
        end
    end

    Rf = hevi_apply_A!(similar(Vf), Vf, params, opfull)
    Rv = hevi_apply_A!(similar(Vv), Vv, params, opvert)

    num = 0.0; den = 0.0
    @inbounds for (q, ieq) in enumerate(opvert.vars)
        qf = opfull.slot[ieq]
        for ip = 1:npoin
            num = max(num, abs(Rf[ip, qf] - Rv[ip, q]))
            den = max(den, abs(Rv[ip, q]))
        end
    end
    rel_vertical = den > 0 ? num / den : num

    # ---- 2. a field varying in x MUST move the horizontal momenta -----------
    # The probe has to be the SAME function of x on every rank, or each rank
    # measures a different wavenumber and the number moves with the partition
    # for reasons that have nothing to do with the operator. extrema(x) is
    # rank-local, so reduce it.
    x  = params.mesh.x
    xmin, xmax = extrema(x)
    let c = get_mpi_comm()
        if MPI.Comm_size(c) > 1
            xmin = MPI.Allreduce(xmin, MPI.MIN, c)
            xmax = MPI.Allreduce(xmax, MPI.MAX, c)
        end
    end
    L  = max(xmax - xmin, eps())
    fill!(Vf, 0.0)
    @inbounds for ip = 1:npoin
        Vf[ip, opfull.slot[5]] = sinpi(2.0 * (x[ip] - xmin) / L)
    end
    Rf2 = hevi_apply_A!(similar(Vf), Vf, params, opfull)
    horiz = 0.0
    @inbounds for ip = 1:npoin
        horiz = max(horiz, abs(Rf2[ip, opfull.slot[2]]))
    end

    comm = get_mpi_comm()
    if MPI.Comm_size(comm) > 1
        rel_vertical = MPI.Allreduce(rel_vertical, MPI.MAX, comm)
        horiz        = MPI.Allreduce(horiz,        MPI.MAX, comm)
    end

    return (rel_vertical = rel_vertical, horiz_response = horiz,
            ok = rel_vertical < 1.0e-10 && horiz > 0.0)
end
