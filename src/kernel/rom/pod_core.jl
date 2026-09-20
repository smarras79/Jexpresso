#---------------------------------------------------------------------------------
# pod_core.jl — Proper Orthogonal Decomposition of a set of snapshots.
#
# THE DECOMPOSITION. Given K snapshots q(x, t_k) of a field on the discrete
# domain, POD is the answer to: which single spatial structure φ(x) captures, on
# average, the most of the signal? Then, which second one does, among those
# orthogonal to the first? And so on. Formally it maximises
#
#       ⟨ |(q', φ)|² ⟩ / (φ, φ)                                            (1)
#
# over φ, where ⟨·⟩ averages over the snapshots and q' = q - q̄ is the
# fluctuation about the temporal mean. The maximiser is the leading eigenfunction
# of the two-point correlation operator, and the whole ordered family {φ_i}
# follows from its spectrum:
#
#       R φ_i = λ_i φ_i ,   R(x,x') = ⟨ q'(x,t) q'(x',t) ⟩ .               (2)
#
# The result is an orthonormal basis ordered by ENERGY: λ_i is the mean square
# of the projection onto φ_i, and no other basis of any given size r captures
# more of ⟨‖q'‖²⟩ than {φ_1 … φ_r}. That optimality is the reason POD is the
# starting point of a reduced-order model: truncating it is the least-damaging
# truncation there is, in the mean-square sense, and the error left behind is
# known exactly in advance,
#
#       ⟨‖q' - Σ_{i≤r} a_i φ_i‖²⟩ / ⟨‖q'‖²⟩ = Σ_{i>r} λ_i / Σ_i λ_i .     (3)
#
# Expanding q(x,t) ≈ q̄(x) + Σ_{i≤r} a_i(t) φ_i(x) and requiring the residual of
# the governing equations to be orthogonal to the retained modes is the Galerkin
# ROM; this file supplies the basis {φ_i}, the coefficients a_i(t_k) and the
# projection/reconstruction operators that such a model is written in terms of.
#
# (Lumley 1967; Sirovich 1987, Q. Appl. Math. 45, 561-571; Holmes, Lumley,
#  Berkooz & Rowley, "Turbulence, Coherent Structures, Dynamical Systems and
#  Symmetry", 2nd ed., CUP 2012.)
#
# THE INNER PRODUCT IS THE MASS MATRIX, and this is not a detail. (1) is posed in
# L²(Ω), so the discrete inner product has to be the discrete L² one,
#
#       (u, v) = ∫_Ω u v dΩ = Σ_ip M_ip u_ip v_ip ,                        (4)
#
# with M the diagonal SEM mass matrix — on the shell, the surface Jacobian times
# the LGL weights, i.e. metrics.M. Using the plain Euclidean product Σ u_ip v_ip
# instead — the default of every off-the-shelf SVD — silently weights every node
# equally, and the LGL nodes are NOT equally spaced: they cluster towards element
# edges like 1/nop², so the edge nodes of every element would count for several
# times their share of the domain, the "modes" would not be orthogonal in L²,
# and (3) would no longer be the error of anything. Here (4) is carried through
# everywhere: in the correlation matrix, in the normalisation of the modes, and
# in the projection that a ROM performs at every time step.
#
# TWO ALGORITHMS, one answer:
#
#   :snapshot  Sirovich's method of snapshots. Form the K×K correlation matrix
#              C_kl = (q'_k, q'_l)/K, take its eigendecomposition, and push the
#              eigenvectors back through the snapshots. Cost O(N K²) with N the
#              number of degrees of freedom, which is what makes POD of a field
#              with N ~ 10⁵ and K ~ 10² a fraction of a second. THE PARALLEL
#              PATH, because C is a sum over nodes: each rank forms its own
#              partial C and one Allreduce of K² numbers completes it. The
#              eigenproblem is then solved redundantly and identically on every
#              rank, and each rank builds its own slice of the modes without any
#              further communication.
#
#   :svd       The singular value decomposition of the mass-weighted snapshot
#              matrix Y = M^{1/2} Q'. Mathematically identical — σ_i² /K = λ_i,
#              and the left singular vectors are M^{1/2} φ_i — but it never forms
#              C and so never squares the condition number. The difference shows
#              up exactly where it matters for a ROM: in the tail of the spectrum,
#              where :snapshot loses the modes below λ_1·ε (ε ≈ 2e-16, so the
#              last ~8 significant digits of the energy) while the SVD resolves
#              down to λ_1·ε². Serial default for that reason.
#
# This file is deliberately free of every Jexpresso type and of MPI: it takes
# arrays and returns a St_pod, and the one thing it needs from a parallel caller
# — a global sum — arrives as a function argument. That is what lets it be
# unit-tested against problems whose POD is known in closed form
# (test/test_pod.jl) without a mesh, a solver or an MPI launcher.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

using LinearAlgebra

export St_pod, pod_from_snapshots
export pod_project, pod_reconstruct, pod_reconstruct!, pod_rank_for_energy,
       pod_truncation_error, pod_mode, pod_inner


"""
    St_pod{TF}

The POD of one field, as returned by [`pod_from_snapshots`](@ref).

| field       | size            | meaning                                          |
|:------------|:----------------|:-------------------------------------------------|
| `name`      |                 | what was decomposed (`"vorticity"`, `"h"`, …)     |
| `comps`     | `ncomp`         | component names, for vector-valued targets        |
| `t`         | `nsnap`         | the snapshot times                                |
| `q̄`         | `npoin × ncomp` | the temporal mean (zero if it was not subtracted) |
| `Φ`         | `npoin × ncomp × nmodes` | the modes, orthonormal in (4)            |
| `λ`         | `nmodes`        | modal energies, descending                        |
| `a`         | `nsnap × nmodes`| temporal coefficients `a_i(t_k)`                  |
| `energy`    | `nmodes`        | `λ_i / Σλ`, the FULL trace in the denominator     |
| `cumenergy` | `nmodes`        | running sum of `energy`                           |
| `total`     |                 | `Σ_i λ_i` over ALL K eigenvalues = `⟨‖q'‖²⟩`      |
| `orthoerr`  |                 | `max |ΦᵀMΦ - I|`, the check that (4) was honoured |

`nmodes ≤ K` because K snapshots span at most a K-dimensional subspace — POD
cannot return more modes than it was shown, whatever the size of the grid.
"""
struct St_pod{TF}
    name      ::String
    comps     ::Vector{String}
    npoin     ::Int
    ncomp     ::Int
    nsnap     ::Int
    t         ::Vector{TF}
    q̄         ::Array{TF,2}
    Φ         ::Array{TF,3}
    λ         ::Vector{TF}
    a         ::Array{TF,2}
    energy    ::Vector{TF}
    cumenergy ::Vector{TF}
    total     ::TF
    lmean     ::Bool
    method    ::Symbol
    orthoerr  ::TF
end

Base.show(io::IO, P::St_pod) = print(io,
    "St_pod(\"", P.name, "\": ", P.npoin, " nodes × ", P.ncomp, " comp, ",
    P.nsnap, " snapshots → ", length(P.λ), " modes, ",
    P.method, ", E₁ = ", round(100*(isempty(P.energy) ? 0.0 : P.energy[1]); digits = 2), "%)")


#
# The two hooks a parallel caller replaces. Serial: nothing to do.
#
_pod_noreduce!(A) = A
_pod_localsign(maxabs, val) = val < 0 ? -one(val) : one(val)


"""
    pod_from_snapshots(X, w, t; kwargs...) -> St_pod

POD of the snapshot set `X`, of size `npoin × ncomp × nsnap`, under the discrete
inner product `(u,v) = Σ_ip w[ip] Σ_c u[ip,c] v[ip,c]` — hand it the diagonal
mass matrix as `w` (see the header, Eq. (4)). `t` carries the snapshot times and
is only stored.

A VECTOR-VALUED target (`ncomp > 1`, e.g. the pair of horizontal velocity
components, or the four conservative variables of the shallow water system) is
decomposed JOINTLY: the components share one set of temporal coefficients, so a
mode is a velocity field rather than two unrelated scalars, and the energy the
spectrum reports is the energy of the vector. That is what a ROM for a coupled
system needs; decomposing the components separately is a different — and for
that purpose wrong — question.

  * `nmodes`   how many modes to keep. `0` (default) keeps every one the
               snapshots can support.
  * `lmean`    subtract the temporal mean first (default `true`). POD of the raw
               field is dominated by a first mode that is essentially the mean,
               which wastes a mode and hides the dynamics; keep it `false` only
               when the mean is genuinely part of what is being modelled.
  * `method`   `:auto` (default), `:svd`, or `:snapshot` — see the header.
  * `rtol`     modes with `λ_i < rtol·λ_1` are discarded as numerical noise.
  * `lparallel`, `reducer!`, `signreduce`
               the parallel hooks: `reducer!(A)` must sum `A` across ranks in
               place, `signreduce(maxabs_local, val_local)` must agree, on every
               rank, on the sign of the global extremum.
"""
function pod_from_snapshots(X::AbstractArray{<:Real,3},
                            w::AbstractVector{<:Real},
                            t::AbstractVector = Float64[];
                            name::String   = "q",
                            comps          = String[],
                            nmodes::Int    = 0,
                            lmean::Bool    = true,
                            method::Symbol = :auto,
                            rtol           = 1.0e-12,
                            lparallel::Bool = false,
                            reducer!        = _pod_noreduce!,
                            signreduce      = _pod_localsign)

    TF = Float64
    npoin, ncomp, nsnap = size(X)

    length(w) >= npoin ||
        error(string(" # ERROR pod_core.jl: the quadrature weights are shorter (", length(w),
                     ") than the snapshots are tall (", npoin, ")."))
    nsnap >= 2 ||
        error(string(" # ERROR pod_core.jl: POD of ", nsnap, " snapshot(s) is not a decomposition.\n",
                     " #   At least two are needed, and a spectrum only means something with more.\n",
                     " #   Raise :ndiagnostics_outputs or :pod_nsnapshots."))

    method in (:auto, :svd, :snapshot) ||
        error(string(" # ERROR pod_core.jl: :pod_method => ", method,
                     " is not one of :auto, :svd, :snapshot."))

    #
    # A rank that owns no node of this field (an empty partition) has w ≡ 0, and
    # M^{1/2} is then not invertible — the SVD path cannot run there. The
    # correlation matrix can: that rank simply contributes nothing to the sum.
    #
    lsvd = if method == :svd
        lparallel && error(" # ERROR pod_core.jl: :pod_method => :svd is serial only; use :snapshot (or :auto) under MPI.")
        true
    elseif method == :snapshot
        false
    else
        !lparallel && all(>(0), @view w[1:npoin])
    end

    #--- the fluctuation matrix Q' = Q - q̄ 1ᵀ, as (npoin·ncomp) × nsnap.
    q̄ = zeros(TF, npoin, ncomp)
    if lmean
        @inbounds for k = 1:nsnap, c = 1:ncomp, ip = 1:npoin
            q̄[ip,c] += TF(X[ip,c,k])
        end
        q̄ ./= nsnap
    end

    Y = Array{TF}(undef, npoin*ncomp, nsnap)
    @inbounds for k = 1:nsnap
        for c = 1:ncomp
            off = (c-1)*npoin
            for ip = 1:npoin
                Y[off+ip, k] = TF(X[ip,c,k]) - q̄[ip,c]
            end
        end
    end

    # the weights, replicated over the components: column-major, so node ip of
    # component c sits at (c-1)*npoin + ip.
    wrep = Vector{TF}(undef, npoin*ncomp)
    @inbounds for c = 1:ncomp, ip = 1:npoin
        wrep[(c-1)*npoin + ip] = TF(w[ip])
    end

    #--- the spectrum, by whichever route
    local λfull::Vector{TF}, Φm::Array{TF,2}, a::Array{TF,2}, total::TF, r::Int

    if lsvd
        sw = sqrt.(wrep)
        F  = svd(sw .* Y)                      # Y ← M^{1/2} Q'
        σ  = F.S
        λfull = (σ .^ 2) ./ nsnap
        total = sum(λfull)
        r     = _pod_rank(λfull, rtol, nmodes)
        Φm    = F.U[:, 1:r] ./ sw              # φ_i = M^{-1/2} u_i
        a     = F.V[:, 1:r] .* transpose(σ[1:r])
    else
        # C = Q'ᵀ M Q' / K, summed over ranks
        C = (transpose(Y) * (wrep .* Y)) ./ nsnap
        C = (C .+ transpose(C)) ./ 2           # exactly symmetric before eigen
        reducer!(C)
        E     = eigen(Symmetric(C))            # ascending
        λfull = reverse(max.(E.values, zero(TF)))
        V     = E.vectors[:, nsnap:-1:1]
        total = sum(λfull)
        r     = _pod_rank(λfull, rtol, nmodes)
        # φ_i = Q' v_i / √(K λ_i)  — the normalisation that makes (φ_i,φ_j) = δ_ij
        scal  = [1.0/sqrt(nsnap*λfull[i]) for i = 1:r]
        Φm    = (Y * V[:, 1:r]) .* transpose(scal)
        a     = V[:, 1:r] .* transpose([sqrt(nsnap*λfull[i]) for i = 1:r])
    end

    #
    # THE SIGN. (2) fixes φ_i only up to ±1 — flipping φ_i and a_i together
    # changes nothing — and LAPACK's choice is not reproducible between the two
    # algorithms, between library versions, or between rank counts. Left alone,
    # a re-run of the same case produces mode plots in inverted colours. Fixing
    # it by "the largest-magnitude entry is positive" is arbitrary but it is a
    # CONVENTION, which is all that is needed.
    #
    @inbounds for i = 1:r
        imax = 1; vmax = zero(TF)
        for n = 1:size(Φm,1)
            if abs(Φm[n,i]) > abs(vmax)
                vmax = Φm[n,i]; imax = n
            end
        end
        s = signreduce(abs(vmax), vmax)
        if s < 0
            @views Φm[:,i] .*= -1
            @views a[:,i]  .*= -1
        end
    end

    #--- orthonormality, ΦᵀMΦ - I, as an honest diagnostic rather than a claim
    G = transpose(Φm) * (wrep .* Φm)
    reducer!(G)
    orthoerr = zero(TF)
    @inbounds for i = 1:r, j = 1:r
        orthoerr = max(orthoerr, abs(G[i,j] - (i == j ? one(TF) : zero(TF))))
    end

    λ         = λfull[1:r]
    energy    = total > 0 ? λ ./ total : zeros(TF, r)
    cumenergy = cumsum(energy)

    Φ = reshape(Φm, npoin, ncomp, r)

    cnames = isempty(comps) ? (ncomp == 1 ? [name] : [string(name, "_", c) for c = 1:ncomp]) :
                              String[string(c) for c in comps]
    length(cnames) == ncomp ||
        error(string(" # ERROR pod_core.jl: ", length(cnames), " component names for ", ncomp, " components."))

    tv = isempty(t) ? collect(TF, 1:nsnap) : TF[TF(tk) for tk in t]
    length(tv) == nsnap ||
        error(string(" # ERROR pod_core.jl: ", length(tv), " times for ", nsnap, " snapshots."))

    return St_pod{TF}(name, cnames, npoin, ncomp, nsnap, tv, q̄, Φ, λ, a,
                      energy, cumenergy, total, lmean,
                      lsvd ? :svd : :snapshot, orthoerr)
end


#
# How many modes the snapshots actually support: everything above the noise
# floor, then capped by what the caller asked for.
#
function _pod_rank(λ::AbstractVector, rtol, nmodes::Int)
    isempty(λ) && return 0
    λ1 = λ[1]
    λ1 > 0 || error(" # ERROR pod_core.jl: every snapshot is identical to the mean — there is nothing to decompose.")
    r = 0
    for i in eachindex(λ)
        λ[i] > rtol*λ1 || break
        r = i
    end
    r = max(r, 1)
    nmodes > 0 && (r = min(r, nmodes))
    return r
end


"""
    pod_inner(P, u, v) -> Real

The discrete L² inner product (4) of two `npoin × ncomp` fields, weighted by the
mass matrix the POD was built with. Needs the weights, which the St_pod does not
store — see [`pod_project`](@ref), which takes them explicitly.
"""
pod_inner(w::AbstractVector, u::AbstractArray, v::AbstractArray) =
    sum(w[ip]*u[ip,c]*v[ip,c] for ip = 1:size(u,1), c = 1:size(u,2))


"""
    pod_project(P, w, q; nmodes = 0, reducer! = identity) -> a

The generalised coordinates of a state `q` (`npoin × ncomp`) in the POD basis,

    a_i = (q - q̄, φ_i) ,   i = 1 … r ,

which is the map a reduced-order model applies to its initial condition, and to
whatever full-order field it has to re-enter. Exact inverse of
[`pod_reconstruct`](@ref) within the span of the retained modes.
"""
function pod_project(P::St_pod{TF}, w::AbstractVector, q::AbstractArray;
                     nmodes::Int = 0, reducer! = _pod_noreduce!) where {TF}
    r = nmodes > 0 ? min(nmodes, length(P.λ)) : length(P.λ)
    size(q,1) == P.npoin && size(q,2) == P.ncomp ||
        error(string(" # ERROR pod_core.jl: pod_project expects a ", P.npoin, " × ", P.ncomp,
                     " state, got ", size(q,1), " × ", size(q,2), "."))
    a = zeros(TF, r)
    @inbounds for i = 1:r, c = 1:P.ncomp, ip = 1:P.npoin
        a[i] += w[ip]*(q[ip,c] - P.q̄[ip,c])*P.Φ[ip,c,i]
    end
    reducer!(a)
    return a
end


"""
    pod_reconstruct(P, a) -> q ;  pod_reconstruct!(q, P, a)

The rank-`r` approximation `q̄ + Σ_{i≤r} a_i φ_i`, i.e. the lift from the reduced
coordinates back to the full grid. `a` may be shorter than the basis, in which
case it truncates.
"""
function pod_reconstruct(P::St_pod{TF}, a::AbstractVector) where {TF}
    q = Array{TF}(undef, P.npoin, P.ncomp)
    return pod_reconstruct!(q, P, a)
end

function pod_reconstruct!(q::AbstractArray, P::St_pod, a::AbstractVector)
    r = min(length(a), length(P.λ))
    @inbounds for c = 1:P.ncomp, ip = 1:P.npoin
        q[ip,c] = P.q̄[ip,c]
    end
    @inbounds for i = 1:r
        ai = a[i]
        for c = 1:P.ncomp, ip = 1:P.npoin
            q[ip,c] += ai*P.Φ[ip,c,i]
        end
    end
    return q
end


"""
    pod_rank_for_energy(P, frac) -> r

The smallest `r` whose modes carry at least `frac` of the total energy — i.e.
how many equations a ROM needs in order to be able to represent that fraction of
the flow. Returns `length(P.λ)` when the retained basis never gets there (the
tail was truncated), so the answer is always a usable rank.
"""
function pod_rank_for_energy(P::St_pod, frac::Real)
    0 < frac <= 1 || error(string(" # ERROR pod_core.jl: an energy fraction of ", frac, " is not in (0,1]."))
    for r in eachindex(P.cumenergy)
        P.cumenergy[r] >= frac && return r
    end
    return length(P.λ)
end


"""
    pod_truncation_error(P) -> e

`e[r+1]` is the relative L² error of the rank-`r` reconstruction, Eq. (3) of the
header, for `r = 0 … nmodes`; `e[1]` is therefore the error of keeping the mean
alone. This is an a-priori curve read off the eigenvalues — no reconstruction is
performed — and it is the lower bound on what ANY r-dimensional linear model can
do with these snapshots.
"""
function pod_truncation_error(P::St_pod{TF}) where {TF}
    n = length(P.λ)
    e = Vector{TF}(undef, n+1)
    P.total > 0 || (fill!(e, zero(TF)); return e)
    tail = P.total
    e[1] = one(TF)
    @inbounds for r = 1:n
        tail -= P.λ[r]
        e[r+1] = sqrt(max(tail, zero(TF))/P.total)
    end
    return e
end


"""
    pod_mode(P, i; comp = 1) -> view

Mode `i`, component `comp`, as a plain `npoin` vector ready to be written to VTK
or handed to the rasterizer.
"""
pod_mode(P::St_pod, i::Int; comp::Int = 1) = view(P.Φ, :, comp, i)
