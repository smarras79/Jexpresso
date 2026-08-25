#=============================================================================
 schur_precond.jl -- the ONE-FIELD column preconditioner for the scalar Schur
                     operator H, built the same way the five-field one is:
                     by PROBING, not by deriving.

 WHY THIS FILE EXISTS
 --------------------
 Step 4 (test/hevi/test_schur_solve.jl) proved the reduction exact: solving the
 scalar system H[P] = rhs and recovering the five fields reproduces the full
 5-field stage solve to 3e-14. The iteration sweep that followed measured H at
 roughly HALF the iterations of the full operator with a column preconditioner
 on both -- but that measurement used, as the preconditioner, the exact inverse
 of the FULL 3D operator restricted to each column's node set. Production does
 not have that. `imex3d_precond!` factorises the VERTICAL-ONLY operator `opv`,
 which is a different matrix: block Jacobi by column keeps the diagonal part of
 the horizontal coupling, and the vertical operator throws all of it away.

 So the honest like-for-like comparison needed a preconditioner that stands in
 the same relation to H as `imex3d_precond!` does to (I - gdt*A). That is what
 this file builds, and the number it produces is the one that should be
 believed over the step-4 sweep.

 THE STENCIL IS TWICE AS WIDE, AND THAT IS THE ONE THING TO GET RIGHT
 --------------------------------------------------------------------
 `assemble_column_band` colours levels modulo S = 2(ngl-1)+1 because A reaches
 ngl-1 levels: the zeta collocation derivative spans an element, DSS links one
 element up and one down. H is not A. It is

     H[P] = P./(beta.*tb) + lam*(W m).grad(tb)./tb + lam*Div[W m],  m = m(P)

 with m = A^-1(-lam*Grad[P] - ...) and A^-1 POINTWISE. Grad reaches ngl-1
 levels and Div reaches ngl-1 more, so H reaches 2(ngl-1) -- and the colouring
 period must be 4(ngl-1)+1, not 2(ngl-1)+1. Using the narrower one does not
 fail loudly: it aliases the outer band entries onto the wrong columns and
 yields a plausible, wrong, still-converging preconditioner. The width is
 therefore ASSERTED against a dense H in test/hevi/test_schur_precond.jl rather
 than argued for here.

 UNLIKE W = I - gdt*A, H IS NOT AN IDENTITY PLUS SOMETHING
 ---------------------------------------------------------
 `assemble_column_band` seeds the band with the identity and subtracts gdt*A,
 because that is literally what W is. H already carries its own zeroth-order
 term (`P./(beta.*tb)`), so this assembly seeds with ZERO and writes the probe
 result straight in. Seeding with the identity here would add a spurious 1 to
 every diagonal entry -- and since the true diagonal is O(1/(beta*tb)) ~ 3e-3
 on the production sounding, that error is ~300x the entry it corrupts and
 would still produce a preconditioner that "works", just a badly wrong one.
=============================================================================#

using LinearAlgebra: LAPACK, BlasInt

"""
    SchurColumnPrecond

One banded LU per owned column of the VERTICAL part of the scalar Schur
operator H, plus the scratch its assembly and application need.

`opv` is a vertical-only (`full = false`) operator carrying ALL FIVE slots and
the advective theta row. All five because `schur_H!` reads the momentum slots
that `schur_grad!` writes and the density slot that `schur_divW!` reads;
advective because the reduction closes on one scalar only with that row -- see
the step-4 header in schur.jl.

`cc` is a `ColumnComm` with `nimp = 1`: this preconditioner moves one number
per level, where the five-field one moves three or five.
"""
struct SchurColumnPrecond{OPT}
    opv::OPT
    st::SchurState{Float64}
    cc::ColumnComm
    topo::ColumnTopology
    n::Int
    kl::Int
    ku::Int
    AB::Vector{Matrix{Float64}}
    ipiv::Vector{Vector{BlasInt}}
    lam::Base.RefValue{Float64}
    nrefactor::Base.RefValue{Int}
    # npoin x 1 staging for column_gather!/scatter!, which take a matrix and a
    # variable list. Kept here so neither the assembly nor the application
    # allocates.
    F::Matrix{Float64}
    P::Vector{Float64}
    HP::Vector{Float64}
end

"""
    schur_column_reach(ngl) -> Int

How many levels above and below a node the scalar Schur operator H couples to:
`2*(ngl-1)`, one `ngl-1` for the gradient and one for the divergence. See the
header, and `test_schur_precond.jl` for the assertion.
"""
schur_column_reach(ngl::Int) = 2 * (ngl - 1)

"""
    assemble_schur_column_band!(pc, params, lam) -> (AB, kl, ku, n)

Probe `schur_H!` and assemble the vertical part of H per owned column, in
LAPACK general-band storage, unfactorised.

`S = 2*reach + 1 = 4*(ngl-1)+1` probes, each of which costs two applications of
the (vertical, hence cheap) operator. At ngl = 5 that is 17 probes against the
9 x nimp the five-field assembly runs -- fewer operator applications than the
five-field preconditioner needs, for a band that is `nimp` times narrower and
`nimp` times shorter.
"""
function assemble_schur_column_band!(pc::SchurColumnPrecond, params, lam::Real)

    npoin = Int(params.mesh.npoin)
    ngl   = Int(params.mesh.ngl)
    topo  = pc.topo
    cc    = pc.cc
    nlev  = topo.nlev
    nown  = cc.nown
    w     = schur_column_reach(ngl)
    S     = 2 * w + 1
    n     = pc.n
    kl    = pc.kl
    ku    = pc.ku
    λ     = Float64(lam)

    lev = level_of_node(topo, npoin)

    diagrow = kl + ku + 1
    AB = pc.AB
    # ZERO, not the identity: H carries its own zeroth-order term. See header.
    @inbounds for ic = 1:nown
        fill!(AB[ic], 0.0)
    end

    P  = pc.P
    HP = pc.HP
    F  = pc.F

    for r = 0:S-1
        fill!(P, 0.0)
        @inbounds for ip = 1:npoin
            l = lev[ip]
            (l != 0 && l % S == r) && (P[ip] = 1.0)
        end

        schur_H!(HP, P, params, pc.opv, λ, pc.st)
        @inbounds for ip = 1:npoin
            F[ip, 1] = HP[ip]
        end
        column_gather!(cc, F, 1:1)

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
                v = cc.X[il, ic]
                v == 0.0 && continue
                Aic[diagrow + il - ilp, ilp] = v
            end
        end
    end

    return AB, kl, ku, n
end

"""
    build_schur_column_precond(params, topo, comm, lam; verify_hook = nothing)
        -> SchurColumnPrecond

Build the vertical-only operator, the one-wide column plan, and the banded LU
of the vertical Schur operator.

`verify_hook(AB, kl, ku, n)` is called with the UNFACTORISED band, for the same
reason `build_column_factorization` offers one: `gbtrf!` overwrites `AB`, so
this is the only moment the assembled matrix can be checked against the
operator it claims to be.
"""
function build_schur_column_precond(params, topo::ColumnTopology, comm::MPI.Comm,
                                    lam::Real; verify_hook = nothing,
                                    lwall_flux::Bool = true)

    npoin = Int(params.mesh.npoin)
    ngl   = Int(params.mesh.ngl)

    # All five slots, vertical only, advective theta row. See the struct docs.
    #
    # `theta_advective` on a `full = false` operator does NOT give it an
    # advective Theta row -- the vertical element kernel has no such branch and
    # keeps the flux one. Its only effect here is to fill `dtbd*`, and that is
    # precisely what is wanted: `schur_H!` reads the CONTINUITY row (through
    # schur_divW!), the MOMENTUM rows (through schur_grad!) and grad(thetabar)
    # pointwise. It never touches the Theta row, so which form that row carries
    # cannot affect the band this file assembles. Verified rather than argued:
    # test_schur_precond.jl matches the assembled band against a dense H_v to
    # 0.000e+00.
    opv = build_hevi_operator(params, topo, [1, 2, 3, 4, 5];
                              lwall_flux = lwall_flux, full = false,
                              theta_advective = true)

    owner, own = assign_column_owners(topo, comm)
    cc = build_column_comm(topo, owner, own, comm, 1)

    w  = schur_column_reach(ngl)
    n  = topo.nlev
    kl = w
    ku = w

    pc = SchurColumnPrecond(opv, SchurState(npoin, opv.nimp), cc, topo,
                            n, kl, ku,
                            [zeros(Float64, 2*kl + ku + 1, n) for _ = 1:cc.nown],
                            [Vector{BlasInt}(undef, n) for _ = 1:cc.nown],
                            Ref(Float64(lam)), Ref(0),
                            zeros(Float64, npoin, 1),
                            zeros(Float64, npoin), zeros(Float64, npoin))

    assemble_schur_column_band!(pc, params, lam)
    verify_hook === nothing || verify_hook(pc.AB, kl, ku, n)
    _factor_schur_band!(pc)

    return pc
end

function _factor_schur_band!(pc::SchurColumnPrecond)
    @inbounds for ic = 1:length(pc.AB)
        _, piv = LAPACK.gbtrf!(pc.kl, pc.ku, pc.n, pc.AB[ic])
        copyto!(pc.ipiv[ic], piv)
    end
    return pc
end

"""
    refactorize_schur!(pc, params, lam)

Rebuild the factorisation for a new `lam`, reusing the storage already there --
the same reason `refactorize!` exists: every diagnostic time is a `tstop`, and
each one moves gamma*dt twice.
"""
function refactorize_schur!(pc::SchurColumnPrecond, params, lam::Real)
    assemble_schur_column_band!(pc, params, lam)
    _factor_schur_band!(pc)
    pc.lam[] = Float64(lam)
    pc.nrefactor[] += 1
    return pc
end

"""
    schur_precond!(P, pc, params, lam) -> P

Apply the inverse of the vertical Schur operator in place to the scalar field
`P`, one banded triangular solve pair per owned column.

Refactorises if `lam` has moved, for the reason in `refactorize_schur!`.
"""
function schur_precond!(P::AbstractVector, pc::SchurColumnPrecond, params, lam::Real)

    # The refactorisation comes FIRST and that is not arbitrary: the assembly
    # probes through this same `pc.F` staging buffer, so staging `P` into it
    # before re-probing would hand the probe the field being preconditioned and
    # get back a band that is not the operator.
    if abs(pc.lam[] - Float64(lam)) > 1.0e-14 * max(abs(Float64(lam)), 1.0)
        _t = hevi_tic()
        refactorize_schur!(pc, params, lam)
        if hevi_prof_on()
            HEVI_PROFILE.t_refac += (time_ns() - _t) * 1e-9; HEVI_PROFILE.n_refac += 1
        end
    end

    # ACCOUNTED TO THE SAME COUNTERS imex3d_precond! uses, so the profile block
    # compares like with like across the two stage solves. Without this the
    # `precond` row reads 0.0 on the Schur path and its cost turns up only as
    # the gap between `sub-accounted` and the solve total -- which is how the
    # first production profile of this path had to be read.
    _tpc = hevi_tic()
    F = pc.F
    @inbounds for ip in eachindex(P)
        F[ip, 1] = P[ip]
    end
    column_gather!(pc.cc, F, 1:1)
    _tb = hevi_tic()
    @inbounds for ic = 1:pc.cc.nown
        LAPACK.gbtrs!('N', pc.kl, pc.ku, pc.n, pc.AB[ic], pc.ipiv[ic],
                      view(pc.cc.X, :, ic))
    end
    if hevi_prof_on()
        HEVI_PROFILE.t_band += (time_ns() - _tb) * 1e-9; HEVI_PROFILE.n_band += 1
    end
    column_scatter!(pc.cc, F, 1:1)
    @inbounds for ip in eachindex(P)
        P[ip] = F[ip, 1]
    end
    if hevi_prof_on()
        HEVI_PROFILE.t_pc += (time_ns() - _tpc) * 1e-9; HEVI_PROFILE.n_pc += 1
    end
    return P
end
