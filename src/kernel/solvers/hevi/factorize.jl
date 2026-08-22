#=============================================================================
 factorize.jl -- one banded LU per vertical column, and the stage solve.

 THE MATRIX IS EXTRACTED, NOT DERIVED
 ------------------------------------
 The operator A applied by hevi_apply_A! couples a node only to nodes within
 +-(ngl-1) levels of it in the same column: inside an element the ζ
 collocation derivative reaches all ngl nodes of the line, and DSS links an
 element only to the ones directly above and below.

 That locality makes the matrix recoverable by *probing*. Colour the levels
 modulo S = 2(ngl-1)+1. A probe that is 1 on one colour and 0 elsewhere hits,
 for each row, exactly one column of the band -- the window of levels a row
 can see is S wide, so it contains exactly one level of any given colour.
 S * nimp applications of A therefore recover every band entry exactly.

 Deriving the same entries by hand would mean re-deriving, and keeping in
 sync, the ζ collocation weights, the metric terms, the DSS sum, the inverse
 mass matrix and the wall flux condition. A discrepancy in any of them would
 not break a convergence test -- the split stays consistent whatever the
 implicit matrix is, because f_exp is computed as rhs! - f_imp. It would just
 make the implicit solve an approximate one, and the vertical acoustic CFL
 would come back with nothing in the output to say why. Probing removes that
 whole failure mode: the matrix is the operator, by construction.

 GAMMA * DT IS CONSTANT, SO THE FACTORISATION IS DONE ONCE
 ---------------------------------------------------------
 Every ARK tableau here has a constant implicit diagonal γ, and the run is
 fixed-step, so W = I - γΔt A never changes. It is factorised once at setup
 and every stage of every step is two triangular solves. If γΔt ever does
 change the factorisation is silently redone and the run says so once.
=============================================================================#

using LinearAlgebra: LAPACK, BlasInt

"""
    ColumnFactorization

Banded LU of `W = I - γΔt A`, one per owned column.

`AB[ic]` is LAPACK general-band storage: `(2kl + ku + 1) x n`, with `W[i,j]`
at `AB[kl + ku + 1 + i - j, j]`, the leading `kl` rows reserved as fill-in
workspace for the factorisation.

`kl = ku = nimp*ngl - 1` because the widest coupling is `ngl-1` levels at
`nimp` variables per level, plus the variables within the level itself.
"""
struct ColumnFactorization
    n::Int
    kl::Int
    ku::Int
    AB::Vector{Matrix{Float64}}
    ipiv::Vector{Vector{BlasInt}}
    gdt::Base.RefValue{Float64}
    nrefactor::Base.RefValue{Int}
end

"""
    level_of_node(topo, npoin) -> Vector{Int}

Level index of every local node; 0 for a node the topology never placed
(which cannot happen on a well-formed mesh, and is worth catching if it does).
"""
function level_of_node(topo::ColumnTopology, npoin::Int)
    lev = zeros(Int, npoin)
    @inbounds for ic = 1:topo.ncol, il = 1:topo.nlev
        ip = topo.node[il, ic]
        ip != 0 && (lev[ip] = il)
    end
    return lev
end

"""
    assemble_column_band(params, op, cc, topo, gdt) -> (AB, kl, ku, n)

Probe the operator and assemble `W = I - γΔt A` per owned column, in LAPACK
general-band storage but not yet factorised -- the setup self-check needs to
multiply by it, which `gbtrf!` would make impossible by overwriting it.
"""
function assemble_column_band(params, op::HEVIOperator, cc::ColumnComm,
                              topo::ColumnTopology, gdt::Real;
                              AB::Union{Nothing, Vector{Matrix{Float64}}} = nothing)

    npoin = Int(params.mesh.npoin)
    ngl   = Int(params.mesh.ngl)
    nimp  = op.nimp
    nlev  = topo.nlev
    nown  = cc.nown
    w     = ngl - 1
    S     = 2 * w + 1
    n     = nimp * nlev
    kl    = nimp * ngl - 1
    ku    = kl

    lev = level_of_node(topo, npoin)

    AB = AB === nothing ? [zeros(Float64, 2*kl + ku + 1, n) for _ = 1:nown] : AB

    # W starts as the identity; the probes subtract γΔt A entry by entry.
    # Zeroing first matters on the reassembly path, where AB is the previous
    # step's LU factors rather than a fresh allocation.
    diagrow = kl + ku + 1
    @inbounds for ic = 1:nown
        fill!(AB[ic], 0.0)
        for j = 1:n
            AB[ic][diagrow, j] = 1.0
        end
    end

    Vprobe = op.V
    Rprobe = op.R

    for r = 0:S-1, m = 1:nimp
        fill!(Vprobe, 0.0)
        @inbounds for ip = 1:npoin
            l = lev[ip]
            (l != 0 && l % S == r) && (Vprobe[ip, m] = 1.0)
        end

        hevi_apply_A!(Rprobe, Vprobe, params, op)
        column_gather!(cc, Rprobe, 1:nimp)

        @inbounds for ic = 1:nown
            Aic = AB[ic]
            for il = 1:nlev
                # the one level of colour r that row `il` can see
                ilp = 0
                for l = max(1, il - w):min(nlev, il + w)
                    if l % S == r
                        ilp = l
                        break
                    end
                end
                ilp == 0 && continue
                j = (ilp - 1) * nimp + m
                for mp = 1:nimp
                    v = cc.X[(il - 1) * nimp + mp, ic]
                    v == 0.0 && continue
                    i = (il - 1) * nimp + mp
                    Aic[diagrow + i - j, j] -= gdt * v
                end
            end
        end
    end

    return AB, kl, ku, n
end

"""
    build_column_factorization(params, op, cc, topo, gdt; verify_hook = nothing)
        -> ColumnFactorization

Assemble the band and factorise it, one banded LU per owned column.

`verify_hook`, when given, is called as `verify_hook(AB, kl, ku, n)` with the
*unfactorised* band, before `gbtrf!` overwrites it. That is the only moment
the assembled matrix can be checked against the operator it is supposed to
be, which is why the check hangs off the assembly rather than living in its
own second pass -- probing is the expensive part of setup and doing it twice
to run a self-check would be a poor trade.
"""
function build_column_factorization(params, op::HEVIOperator, cc::ColumnComm,
                                    topo::ColumnTopology, gdt::Real;
                                    verify_hook = nothing)

    AB, kl, ku, n = assemble_column_band(params, op, cc, topo, gdt)
    verify_hook === nothing || verify_hook(AB, kl, ku, n)

    ipiv = [Vector{BlasInt}(undef, n) for _ = 1:length(AB)]
    @inbounds for ic = 1:length(AB)
        _, piv = LAPACK.gbtrf!(kl, ku, n, AB[ic])
        copyto!(ipiv[ic], piv)
    end

    return ColumnFactorization(n, kl, ku, AB, ipiv, Ref(Float64(gdt)), Ref(0))
end

"""
    refactorize!(fac, params, op, cc, topo, gdt)

Rebuild the factorisation for a new γΔt, reusing the storage already there.

This is NOT an exceptional path. OrdinaryDiffEq shortens the step that lands
on a `tstop`, and Jexpresso passes every diagnostic and LES-statistics time as
a tstop -- on the LESICP2 deck that is ~180 of them, each costing two changes
of γΔt (one onto the short step, one back). Allocating a fresh set of band
matrices each time would churn tens of megabytes per event for no reason, so
the assembly writes into the existing arrays.

The cost itself is small: S*nimp applications of the (cheap) operator plus one
banded LU per column, together on the order of a single time step. At ~360
events in a 200 000-step run that is around 0.2% overhead.
"""
function refactorize!(fac::ColumnFactorization, params, op::HEVIOperator,
                      cc::ColumnComm, topo::ColumnTopology, gdt::Real)

    assemble_column_band(params, op, cc, topo, gdt; AB = fac.AB)
    @inbounds for ic = 1:length(fac.AB)
        _, piv = LAPACK.gbtrf!(fac.kl, fac.ku, fac.n, fac.AB[ic])
        copyto!(fac.ipiv[ic], piv)
    end
    fac.gdt[] = Float64(gdt)
    fac.nrefactor[] += 1
    return fac
end

"""
    hevi_column_solve!(dst, src, params, gdt)

The stage solve: `dst = qref + (I - γΔt A)^-1 (src - qref)`, with `dst == src`
in every equation that is not implicit.

`src` and `dst` are the flat `npoin*neqs` state vectors the rest of the code
uses. Equations outside `op.vars` are copied straight through, which is
exactly right: `A` has zero rows there, so the stage equation for them is the
identity.
"""
function hevi_column_solve!(dst, src, params, gdt::Real)

    hevi = params.hevi
    op   = hevi.op
    cc   = hevi.cc
    fac  = hevi.fac
    npoin = Int(params.mesh.npoin)
    neqs  = Int(params.neqs)
    qe    = params.qp.qe

    # γΔt changes whenever the integrator shortens a step to land on a tstop,
    # so this fires a couple of hundred times in a normal run -- see
    # refactorize!. Solving with a stale W instead would be cheaper and wrong,
    # and would present as an unexplained stability limit rather than as a
    # mistake, which is exactly the failure mode worth paying to avoid.
    if abs(fac.gdt[] - Float64(gdt)) > 1.0e-14 * max(abs(gdt), 1.0)
        hevi_trace("    refactorising: gdt ", fac.gdt[], " -> ", Float64(gdt),
                   " (rebuild #", fac.nrefactor[] + 1, ")")
        _t = hevi_tic()
        refactorize!(fac, params, op, cc, hevi.topo, gdt)
        if hevi_prof_on()
            HEVI_PROFILE.t_refac += (time_ns() - _t) * 1e-9; HEVI_PROFILE.n_refac += 1
        end
        hevi_trace("    refactorisation done")
    end

    # deviation from the reference state, packed to the implicit variables
    V = op.V
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            V[ip, q] = src[off + ip] - qe[ip, ieq]
        end
    end

    hevi_trace("    gather")
    column_gather!(cc, V, 1:op.nimp)

    hevi_trace("    ", cc.nown, " banded solves")
    @inbounds for ic = 1:cc.nown
        LAPACK.gbtrs!('N', fac.kl, fac.ku, fac.n, fac.AB[ic], fac.ipiv[ic],
                      view(cc.X, :, ic))
    end

    hevi_trace("    scatter")
    column_scatter!(cc, V, 1:op.nimp)

    dst === src || copyto!(dst, src)
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            dst[off + ip] = V[ip, q] + qe[ip, ieq]
        end
    end
    return dst
end
