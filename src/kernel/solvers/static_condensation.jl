#=============================================================================
 static_condensation.jl — the static-condensation SOLVE of a symmetric SEM
 system (elementLearning_Axb! with :lstatic_condensation, 2D and 3D), as a
 threaded kernel that keeps nothing per element:

   B_uu u_u = f_u − A_uΓ g − Σ_e A_bo A_oo⁻¹ (f_o − A_oΓ g)        skeleton
   B_uu     = A_uu − Σ_e A_bo A_oo⁻¹ A_ob                         Schur complement
   u_o      = A_oo⁻¹ (f_o − A_ob u_b)                             per element

 (u: skeleton unknowns = ∂O, less the pinned first one when there is no
 Dirichlet set Γ; o: an element's interior nodes; b: its boundary nodes.)

 Per element, with A_oo = L Lᵀ (LAPACK potrf) and W = L⁻¹ A_ob (trsm):
   A_bo A_oo⁻¹ A_ob = Wᵀ W   (syrk),   A_bo A_oo⁻¹ v = Wᵀ (L⁻¹ v)   (trsv, gemv)
 against the legacy loop's dense inverse (≈2 n³ instead of n³/3) and its
 scalar triple loop for the product.

 B is assembled straight into its final sparsity pattern (built once,
 symbolically, from the elements sharing each skeleton node), not from
 per-element triplets + sparse(): at 16³ elements, N = 8, those triplets
 were 6·10⁸ (14 GB) before sorting. The elements are coloured so that no two
 of one colour share a skeleton node; each colour is condensed in parallel,
 conflict-free and in a fixed order, so the result is deterministic and B is
 EXACTLY symmetric (entry (r,c) and (c,r) receive the same values in the same
 order): no symmetrisation copy, no pinned-submatrix copy (the pinned node is
 left out of the unknowns, i.e. treated as Dirichlet with g = 0).

 Nothing is stored per element: the interior recovery re-extracts A_oo, A_ob
 and refactors A_oo (n³/3, a tenth of the condensation's work) instead of
 keeping L and A_ob for every element (2·nelem·n_o·(n_o+n_b) values: 8 GB at
 16³ elements, N = 8).

 Applicability (checked; otherwise the legacy loop runs): every element's
 boundary nodes are skeleton or Dirichlet nodes, no node repeats within an
 element, its interior nodes belong to it alone, and A is exactly symmetric
 on the non-Dirichlet rows and columns.
=============================================================================#

# Per-thread work arrays of one element (sizes n_o = interior, n_b = boundary)
struct _SCWork
    loc  :: Vector{Int32}      # global node -> local index in the current element (0 elsewhere)
    Aoo  :: Matrix{Float64}    # A_oo, then its Cholesky factor L (lower)
    W    :: Matrix{Float64}    # A_ob, then L⁻¹ A_ob
    S    :: Matrix{Float64}    # Wᵀ W (upper triangle)
    v    :: Vector{Float64}    # f_o − A_ob g_b, then L⁻¹ (…)
    gb   :: Vector{Float64}    # boundary values of the element
    dv   :: Vector{Float64}    # Wᵀ L⁻¹ (…)
    key  :: Vector{Int}        # (skeleton index << 16) | active index, sorted
    act  :: Vector{Int}        # the element's active boundary nodes (local indices)
    amap :: Vector{Int}        # local boundary index -> active index (0: inactive)
end
_SCWork(n, no, nb) = _SCWork(zeros(Int32, n), zeros(no, no), zeros(no, nb), zeros(nb, nb),
                             zeros(no), zeros(nb), zeros(nb), Vector{Int}(undef, nb),
                             Vector{Int}(undef, nb), zeros(Int, nb))

# Run f(chunk, lo, hi) over 1:n in nchunks contiguous chunks, in parallel
@inline function _sc_chunks(f::F, n::Int, nchunks::Int) where {F}
    Threads.@threads :static for c = 1:nchunks
        lo = div((c - 1) * n, nchunks) + 1
        hi = div(c * n, nchunks)
        f(c, lo, hi)
    end
    return nothing
end

# Does column `col` of A hold row `row`? -> its position in A.nzval, or 0
@inline function _sc_findrow(A::SparseMatrixCSC, row::Integer, col::Integer)
    r1 = Int(A.colptr[col]); r2 = Int(A.colptr[col+1]) - 1
    p  = searchsortedfirst(A.rowval, row, r1, r2, Base.Order.Forward)
    return p <= r2 && A.rowval[p] == row ? p : 0
end

"""
    el_sc_solve!(u, A, f, conn, nb, ∂O, Γ, gΓ; solver = :direct, amg_method = "sa",
                 amg_rtol = 1e-12, ordering = :cholmod) -> Bool

Static-condensation solve of `A u = f` (see the file header). `conn` is the
element connectivity with each element's `nb` boundary nodes first; `∂O` the
internal skeleton, `Γ` the Dirichlet nodes with values `gΓ` (empty: the
system is singular with the constants in its null space and ∂O[1] is pinned
to 0). Writes the solution into `u` and returns true; returns false, having
done nothing, if the kernel does not apply (then use the legacy loop).
Records :sc_extract (maps, checks, pattern of B, colouring), :sc_condense,
:sc_factor, :sc_skeleton and :sc_recover in JX_TIMINGS.
"""
function el_sc_solve!(u::AbstractArray, A::SparseMatrixCSC{Float64}, f::AbstractVector,
                      conn::AbstractMatrix{<:Integer}, nb::Int,
                      ∂O::AbstractVector{<:Integer}, Γ::AbstractVector{<:Integer},
                      gΓ::AbstractVector;
                      solver = :direct, amg_method = "sa", amg_rtol::Real = 1e-12,
                      ordering::Symbol = :cholmod)
    t0 = time_ns()
    nt = Threads.nthreads()
    R  = _sc_roles(A, conn, nb, ∂O, Γ, gΓ, nt)
    R === nothing && return false
    (; n, nelem, npel, no, unk, ns, pos, own, ug, sgn) = R

    # ── ACTIVE boundary nodes: those coupled (in A) to the element's interior.
    #    Only they enter A_bo A_oo⁻¹ A_ob: with collocation (Inexact) the
    #    interior couples only to the face-interior nodes, not to the element's
    #    edges and vertices, whose pairs would be stored zeros in B (an AMG
    #    hierarchy and CG matvecs twice as expensive at N = 5).
    active = Matrix{Bool}(undef, nb, nelem)
    _sc_chunks(nelem, nt) do c, lo, hi
        @inbounds for e = lo:hi, j = 1:nb
            a = false
            for p in nzrange(A, conn[e, j])
                own[A.rowval[p]] == e && (a = true;  break)
            end
            active[j, e] = a
        end
        return nothing
    end
    # ── skeleton unknown -> the elements on which it is active ───────────────
    eptr, elist = _sc_element_lists(conn, pos, active, nb, ns)
    B      = _sc_pattern(A, conn, pos, active, unk, eptr, elist, nb, ns, nt)
    colors = _sc_colour(conn, pos, active, eptr, elist, nb, nelem)
    JX_TIMINGS[:sc_extract] = (time_ns() - t0) / 1e9

    # ── condensation ──────────────────────────────────────────────────────────
    t0  = time_ns()
    rhs = zeros(Float64, ns)
    @inbounds for k = 1:ns;  rhs[k] = sgn * f[unk[k]];  end
    _sc_add_A!(B, rhs, A, unk, pos, Γ, gΓ, sgn, nt)
    work = [_SCWork(n, no, nb) for _ = 1:nt]
    nblas = BLAS.get_num_threads()
    nt > 1 && BLAS.set_num_threads(1)               # one BLAS thread per Julia thread
    try
        for (cp1, cp2) in colors
            ne = cp2 - cp1 + 1
            _sc_chunks(ne, min(nt, ne)) do c, lo, hi
                w = work[c]
                for q = cp1+lo-1:cp1+hi-1
                    _sc_condense_element!(B, rhs, w, A, f, ug, conn, pos, active, Int(colors.elems[q]), nb, no, sgn)
                end
            end
        end
    finally
        BLAS.set_num_threads(nblas)
    end
    JX_TIMINGS[:sc_condense] = (time_ns() - t0) / 1e9

    # ── skeleton solve ────────────────────────────────────────────────────────
    us = _el_skeleton_solve_spd(B, rhs; solver = solver, amg_method = amg_method,
                                amg_rtol = amg_rtol, ordering = ordering)

    # ── interior recovery ─────────────────────────────────────────────────────
    t0 = time_ns()
    @inbounds for k = 1:ns;  ug[unk[k]] = us[k];  end
    nt > 1 && BLAS.set_num_threads(1)
    try
        _sc_chunks(nelem, min(nt, nelem)) do c, lo, hi
            w = work[c]
            for e = lo:hi
                _sc_recover_element!(ug, w, A, f, conn, active, e, nb, no, sgn)
            end
        end
    finally
        BLAS.set_num_threads(nblas)
    end
    @inbounds for i = 1:n;  u[i] = ug[i];  end
    JX_TIMINGS[:sc_recover] = (time_ns() - t0) / 1e9
    return true
end

# Roles of the nodes and the checks of el_sc_solve! / el_sc_schur_cg!:
#   unk  skeleton unknowns (∂O, less the pinned ∂O[1] when Γ is empty)
#   pos  node -> index in unk (0 elsewhere);  own  interior node -> its element
#   ug   the solution vector, holding the fixed values (Γ: gΓ; pinned: 0)
#   sgn  ±1 such that sgn·A_oo is positive definite (the SEM "L" may be −K)
# `nothing` if the kernels do not apply.
function _sc_roles(A::SparseMatrixCSC, conn, nb, ∂O, Γ, gΓ, nt)
    n  = size(A, 1)
    nelem, npel = size(conn)
    no = npel - nb
    (no > 0 && nb > 0 && length(∂O) > 1) || return nothing
    unk  = isempty(Γ) ? ∂O[2:end] : collect(∂O)
    ns   = length(unk)
    pos  = zeros(Int32, n)
    role = zeros(Int8, n)                    # 1 unknown, 2 fixed, 3 interior
    own  = zeros(Int32, n)
    ug   = zeros(Float64, n)
    @inbounds for (k, g) in enumerate(unk);  pos[g] = k;  role[g] = 1;  end
    @inbounds for g in ∂O;  role[g] == 0 && (role[g] = 2);  end   # the pinned node
    @inbounds for (i, g) in enumerate(Γ);  role[g] = 2;  ug[g] = gΓ[i];  end
    @inbounds for e = 1:nelem
        for j = 1:nb
            role[conn[e, j]] in (1, 2) || return nothing   # boundary node off the skeleton
        end
        for j = nb+1:npel
            g = conn[e, j]
            role[g] == 0 || return nothing                  # interior node shared or on the skeleton
            role[g] = 3;  own[g] = e
        end
    end
    any(iszero, role) && return nothing
    stamp = zeros(Int32, n)                                 # no node twice in an element
    @inbounds for e = 1:nelem, j = 1:npel
        g = conn[e, j];  stamp[g] == e && return nothing;  stamp[g] = e
    end
    _sc_symmetric(A, role, nt) || return nothing
    g1 = conn[1, nb+1];  p1 = _sc_findrow(A, g1, g1)
    sgn = p1 > 0 && A.nzval[p1] < 0 ? -1.0 : 1.0
    return (; n, nelem, npel, no, unk, ns, pos, role, own, ug, sgn)
end

# A exactly symmetric on the rows and columns that are not fixed (role 2)?
function _sc_symmetric(A::SparseMatrixCSC, role::Vector{Int8}, nt::Int)
    n  = size(A, 2)
    ok = fill(true, nt)
    _sc_chunks(n, nt) do c, lo, hi
        @inbounds for j = lo:hi
            role[j] == 2 && continue
            for p in nzrange(A, j)
                i = A.rowval[p]
                role[i] == 2 && continue
                q = _sc_findrow(A, j, i)
                if (q == 0 ? 0.0 : A.nzval[q]) != A.nzval[p]
                    ok[c] = false;  return nothing
                end
            end
        end
        return nothing
    end
    return all(ok)
end

# CSR-like lists: for skeleton unknown k, the elements elist[eptr[k]:eptr[k+1]-1]
# on which it is active
function _sc_element_lists(conn, pos, active, nb, ns)
    nelem = size(conn, 1)
    eptr  = zeros(Int, ns + 1)
    @inbounds for e = 1:nelem, j = 1:nb
        k = pos[conn[e, j]];  k > 0 && active[j, e] && (eptr[k+1] += 1)
    end
    eptr[1] = 1
    @inbounds for k = 1:ns;  eptr[k+1] += eptr[k];  end
    elist = Vector{Int32}(undef, eptr[ns+1] - 1)
    next  = eptr[1:ns]
    @inbounds for e = 1:nelem, j = 1:nb
        k = pos[conn[e, j]]
        k > 0 && active[j, e] || continue
        elist[next[k]] = e;  next[k] += 1
    end
    return eptr, elist
end

# Symbolic pattern of B (ns × ns): column k holds every skeleton unknown of
# every element on k (plus A's own entries); rows sorted, values zero.
function _sc_pattern(A, conn, pos, active, unk, eptr, elist, nb, ns, nt)
    colptr = zeros(Int, ns + 1)
    marks  = [zeros(Int32, ns) for _ = 1:nt]
    _sc_chunks(ns, nt) do c, lo, hi                   # pass 1: count
        for k = lo:hi
            colptr[k+1] = _sc_column!(nothing, 0, marks[c], k, A, conn, pos, active, unk, eptr, elist, nb)
        end
        return nothing
    end
    colptr[1] = 1
    @inbounds for k = 1:ns;  colptr[k+1] += colptr[k];  end
    for mk in marks;  fill!(mk, Int32(0));  end
    rowval = Vector{Int}(undef, colptr[ns+1] - 1)
    _sc_chunks(ns, nt) do c, lo, hi                   # pass 2: write, sort
        for k = lo:hi
            p0 = colptr[k] - 1
            m  = _sc_column!(rowval, p0, marks[c], k, A, conn, pos, active, unk, eptr, elist, nb)
            sort!(view(rowval, p0+1:p0+m); alg = Base.Sort.QuickSort)
        end
        return nothing
    end
    return SparseMatrixCSC(ns, ns, colptr, rowval, zeros(Float64, length(rowval)))
end

# The distinct rows of column k of B (written at rowval[p0+1:…] unless
# rowval === nothing); returns their number
@inline function _sc_column!(rowval, p0, mark, k, A, conn, pos, active, unk, eptr, elist, nb)
    kk = Int32(k);  m = 0
    @inbounds for q = eptr[k]:eptr[k+1]-1
        e = elist[q]
        for j = 1:nb
            active[j, e] || continue
            r = pos[conn[e, j]]
            (r == 0 || mark[r] == kk) && continue
            mark[r] = kk;  m += 1
            rowval === nothing || (rowval[p0+m] = r)
        end
    end
    @inbounds for p in nzrange(A, unk[k])
        r = pos[A.rowval[p]]
        (r == 0 || mark[r] == kk) && continue
        mark[r] = kk;  m += 1
        rowval === nothing || (rowval[p0+m] = r)
    end
    return m
end

# Greedy colouring: no two elements of one colour share an active skeleton unknown.
# Returns the elements grouped by colour, iterable as (first, last) ranges.
struct _SCColours
    elems :: Vector{Int32}
    ptr   :: Vector{Int}
end
Base.length(c::_SCColours) = length(c.ptr) - 1
Base.iterate(c::_SCColours, i = 1) = i > length(c) ? nothing : ((c.ptr[i], c.ptr[i+1] - 1), i + 1)

function _sc_colour(conn, pos, active, eptr, elist, nb, nelem)
    col  = zeros(Int32, nelem)
    seen = zeros(Int, nelem + 1)                 # seen[c] == e: colour c is taken next to e
    ncol = 0
    @inbounds for e = 1:nelem
        for j = 1:nb
            k = pos[conn[e, j]];  (k == 0 || !active[j, e]) && continue
            for q = eptr[k]:eptr[k+1]-1
                c2 = col[elist[q]];  c2 > 0 && (seen[c2] = e)
            end
        end
        c = 1
        while seen[c] == e;  c += 1;  end
        col[e] = c;  ncol = max(ncol, c)
    end
    ptr = zeros(Int, ncol + 1);  ptr[1] = 1
    @inbounds for e = 1:nelem;  ptr[col[e]+1] += 1;  end
    for c = 1:ncol;  ptr[c+1] += ptr[c];  end
    next  = ptr[1:ncol]
    elems = Vector{Int32}(undef, nelem)
    @inbounds for e = 1:nelem
        elems[next[col[e]]] = e;  next[col[e]] += 1
    end
    return _SCColours(elems, ptr)
end

# B ← sgn·A_uu (its pattern holds A's), rhs −= sgn·A_uΓ g; column-owned: no conflicts
function _sc_add_A!(B, rhs, A, unk, pos, Γ, gΓ, sgn, nt)
    ns = length(unk)
    _sc_chunks(ns, nt) do c, lo, hi
        @inbounds for k = lo:hi
            g = unk[k]
            for p in nzrange(A, g)
                r = pos[A.rowval[p]];  r == 0 && continue
                B.nzval[_sc_findrow(B, r, k)] += sgn * A.nzval[p]
            end
        end
        return nothing
    end
    @inbounds for (i, γ) in enumerate(Γ)
        gv = gΓ[i];  iszero(gv) && continue
        for p in nzrange(A, γ)
            r = pos[A.rowval[p]];  r == 0 && continue
            rhs[r] -= sgn * A.nzval[p] * gv
        end
    end
    return nothing
end

# A_oo -> w.Aoo and the active columns of A_ob -> w.W[:, 1:na] (scaled by
# sgn), from the element's columns of A; returns na
@inline function _sc_extract!(w::_SCWork, A, conn, active, e, nb, no, sgn)
    npel = nb + no
    na = 0
    @inbounds for j = 1:nb
        if active[j, e]
            na += 1;  w.act[na] = j;  w.amap[j] = na
        else
            w.amap[j] = 0
        end
    end
    @inbounds for j = 1:npel;  w.loc[conn[e, j]] = j;  end
    fill!(w.Aoo, 0.0);  fill!(view(w.W, :, 1:na), 0.0)
    @inbounds for lc = 1:npel
        ac = lc > nb ? 0 : w.amap[lc]
        lc <= nb && ac == 0 && continue
        for p in nzrange(A, conn[e, lc])
            lr = Int(w.loc[A.rowval[p]]) - nb
            lr > 0 || continue                       # interior rows only
            if lc > nb
                w.Aoo[lr, lc-nb] = sgn * A.nzval[p]
            else
                w.W[lr, ac] = sgn * A.nzval[p]
            end
        end
    end
    @inbounds for j = 1:npel;  w.loc[conn[e, j]] = 0;  end
    return na
end

@inline function _sc_potrf!(Aoo)
    _, info = LAPACK.potrf!('L', Aoo)
    info == 0 || error(" # static condensation: an element's interior block A_oo is not " *
                       "positive definite (LAPACK potrf info = $info).")
    return nothing
end

# v = sgn f_o − A_ob g_b over the active boundary nodes; returns v
@inline function _sc_load!(w::_SCWork, f, ug, conn, e, nb, no, na, sgn, skipzero::Bool)
    nz = false
    @inbounds for a = 1:na
        w.gb[a] = ug[conn[e, w.act[a]]];  nz |= !iszero(w.gb[a])
    end
    @inbounds for i = 1:no;  w.v[i] = sgn * f[conn[e, nb+i]];  end
    (nz || !skipzero) && na > 0 &&
        BLAS.gemv!('N', -1.0, view(w.W, :, 1:na), view(w.gb, 1:na), 1.0, w.v)
    return w.v
end

# One element's Schur complement and load correction, scattered into B and rhs
function _sc_condense_element!(B, rhs, w::_SCWork, A, f, ug, conn, pos, active, e, nb, no, sgn)
    na = _sc_extract!(w, A, conn, active, e, nb, no, sgn)
    _sc_load!(w, f, ug, conn, e, nb, no, na, sgn, true)   # f_o − A_ob g_b (g: Dirichlet values)
    _sc_potrf!(w.Aoo)                                     # A_oo = L Lᵀ
    na == 0 && return nothing
    W = view(w.W, :, 1:na);  S = view(w.S, 1:na, 1:na);  dv = view(w.dv, 1:na)
    BLAS.trsm!('L', 'L', 'N', 'N', 1.0, w.Aoo, W)         # W = L⁻¹ A_ob
    BLAS.trsv!('L', 'N', 'N', w.Aoo, w.v)                 # v = L⁻¹ v
    BLAS.syrk!('U', 'T', 1.0, W, 0.0, S)                  # S = Wᵀ W = A_bo A_oo⁻¹ A_ob
    BLAS.gemv!('T', 1.0, W, w.v, 0.0, dv)                 # dv = A_bo A_oo⁻¹ v

    # the element's active skeleton unknowns, sorted by their index in B
    m = 0
    @inbounds for a = 1:na
        k = Int(pos[conn[e, w.act[a]]]);  k == 0 && continue
        m += 1;  w.key[m] = (k << 16) | a
    end
    key = view(w.key, 1:m)
    sort!(key; alg = Base.Sort.QuickSort)
    @inbounds for c = 1:m
        k = key[c] >> 16;  j = key[c] & 0xffff
        rhs[k] -= dv[j]
        # column k of B: walk its sorted rows along the element's sorted unknowns
        p = B.colptr[k]
        for b = 1:m
            r = key[b] >> 16;  i = key[b] & 0xffff
            while B.rowval[p] < r;  p += 1;  end
            B.nzval[p] -= i <= j ? S[i, j] : S[j, i]
        end
    end
    return nothing
end

# One element's interior values u_o = A_oo⁻¹ (f_o − A_ob u_b)
function _sc_recover_element!(ug, w::_SCWork, A, f, conn, active, e, nb, no, sgn)
    na = _sc_extract!(w, A, conn, active, e, nb, no, sgn)
    _sc_load!(w, f, ug, conn, e, nb, no, na, sgn, false)
    _sc_potrf!(w.Aoo)
    BLAS.trsv!('L', 'N', 'N', w.Aoo, w.v)
    BLAS.trsv!('L', 'T', 'N', w.Aoo, w.v)
    @inbounds for i = 1:no;  ug[conn[e, nb+i]] = w.v[i];  end
    return nothing
end

"""
    _el_skeleton_solve_spd(B, rhs; solver, amg_method, amg_rtol, ordering) -> u

Skeleton solve for an exactly symmetric positive definite B (el_sc_solve!):
sparse Cholesky (CHOLMOD; `ordering = :metis` for METIS nested dissection)
or AMG-preconditioned CG, with no symmetrisation or pinning copies.
Records :sc_factor and :sc_skeleton.
"""
function _el_skeleton_solve_spd(B::SparseMatrixCSC, rhs::Vector{Float64};
                                solver = :direct, amg_method = "sa", amg_rtol = 1e-12,
                                ordering::Symbol = :cholmod)
    s = Symbol(lowercase(string(solver)))
    s in (:direct, :amg) ||
        error(" # el_skeleton_solve: :EL_skeleton_solver => \"$solver\"; expected \"direct\" or \"amg\".")
    if s === :direct
        F = jx_phase(:sc_factor) do
                _el_skeleton_cholesky(B; perm = ordering === :metis ? jx_metis_perm(B) : nothing)
            end
        return jx_phase(() -> F \ rhs, :sc_skeleton)
    else
        S = jx_phase(() -> jx_amg_setup(B; method = amg_method), :sc_factor)
        return jx_phase(() -> jx_amg_solve(S, rhs; rtol = amg_rtol), :sc_skeleton)
    end
end


#=============================================================================
 Tensor-product, matrix-free static condensation for AMG-CG (3D, affine
 hexahedra): el_sc_schur_cg!

 On an affine (box) element of the collocation SEM the interior block is a
 sum of Kronecker products of the 1D GLL stiffness K̂ and weights ω̂:
     sgn·A_oo = c_x  ω̂⊗ω̂⊗K̂  +  c_y  ω̂⊗K̂⊗ω̂  +  c_z  K̂⊗ω̂⊗ω̂      (i fastest)
 (interior rows/columns 2:N of K̂, ω̂). With the 1D generalized eigenpairs
 K̂ s = λ ω̂ s, Sᵀ ω̂ S = I, FAST DIAGONALIZATION inverts it in O(N⁴):
     A_oo⁻¹ = (S⊗S⊗S) diag(1/(c_x λ_i + c_y λ_j + c_z λ_k)) (S⊗S⊗S)ᵀ ,
 with three numbers per element (c_x, c_y, c_z), read off A and checked
 against every entry of A_oo (any other element: the kernel does not apply).

 The Schur complement is then never formed. CG runs on the skeleton with
     S p = ( Â [p; −Â_oo⁻¹ Â_ob p] )_b            (one sparse matvec + FDM)
 preconditioned by the AMG V-cycle of the FULL system Â (that of the SEM
 AMG-CG) restricted to the skeleton, (P⁻¹)_bb: since (Â⁻¹)_bb = S⁻¹ exactly,
 it is spectrally equivalent to S⁻¹ with P's constants. No B (10× the
 nonzeros of A at N = 8), no AMG hierarchy of a dense skeleton, nothing
 stored per element but (c_x, c_y, c_z).
=============================================================================#

# FDM data of the interior (m = N−1 points per direction)
struct _SCFDM
    m  :: Int
    S  :: Matrix{Float64}       # eigenvectors, Sᵀ ω̂ S = I
    St :: Matrix{Float64}       # Sᵀ
    λ  :: Vector{Float64}
    c  :: Matrix{Float64}       # (3, nelem): c_x, c_y, c_z
end

# Per-thread FDM / matvec work arrays
struct _SCFDMWork
    x  :: Vector{Float64}       # m³
    t  :: Vector{Float64}       # m³
    Aoo :: Matrix{Float64}      # m³ × m³, for the check of c only
    loc :: Vector{Int32}
end

# 1D GLL stiffness K̂[a,b] = Σ_k ω_k ψ'_a(ξ_k) ψ'_b(ξ_k)  (dψ[a,k] = ψ'_a(ξ_k))
sc_gll_stiffness(dψ::AbstractMatrix, ω::AbstractVector) =
    [sum(ω[k] * dψ[a, k] * dψ[b, k] for k in eachindex(ω)) for a in axes(dψ, 1), b in axes(dψ, 1)]

# X ← A_oo⁻¹ X on one element (X, t: m³ vectors, i fastest)
@inline function _sc_fdm_solve!(X::Vector{Float64}, t::Vector{Float64}, F::_SCFDM, e::Int)
    m = F.m;  m2 = m * m
    _sc_modes!(t, X, F.St, m)                                # t = (Sᵀ⊗Sᵀ⊗Sᵀ) X
    cx = F.c[1, e];  cy = F.c[2, e];  cz = F.c[3, e];  λ = F.λ
    @inbounds for k = 1:m, j = 1:m
        d0 = cy * λ[j] + cz * λ[k];  o = (k - 1) * m2 + (j - 1) * m
        @simd for i = 1:m
            t[o+i] /= cx * λ[i] + d0
        end
    end
    _sc_modes!(X, t, F.S, m)                                 # X = (S⊗S⊗S) t
    return X
end

# Y ← (Q⊗Q⊗Q) X for an m×m×m tensor (i fastest), in place of Y; X is used as scratch
@inline function _sc_modes!(Y::Vector{Float64}, X::Vector{Float64}, Q::Matrix{Float64}, m::Int)
    m2 = m * m
    # direction i: Y[:, jk] = Q X[:, jk]
    @inbounds for jk = 0:m2-1
        o = jk * m
        for a = 1:m
            s = 0.0
            @simd for b = 1:m;  s += Q[a, b] * X[o+b];  end
            Y[o+a] = s
        end
    end
    # direction j: X[i, :, k] = Q Y[i, :, k]
    @inbounds for k = 0:m-1, a = 1:m
        ok = k * m2
        for i = 1:m;  X[ok+(a-1)*m+i] = 0.0;  end
        for b = 1:m
            q = Q[a, b];  ob = ok + (b - 1) * m;  oa = ok + (a - 1) * m
            @simd for i = 1:m;  X[oa+i] += q * Y[ob+i];  end
        end
    end
    # direction k: Y[ij, :] = Q X[ij, :]
    @inbounds for a = 1:m
        oa = (a - 1) * m2
        for ij = 1:m2;  Y[oa+ij] = 0.0;  end
        for b = 1:m
            q = Q[a, b];  ob = (b - 1) * m2
            @simd for ij = 1:m2;  Y[oa+ij] += q * X[ob+ij];  end
        end
    end
    return Y
end

# (c_x, c_y, c_z) of every element from A, each A_oo checked entry by entry
# against the Kronecker form; nothing if an element is not of that form
function _sc_fdm_setup(A, conn, nb, no, sgn, K̂, ω, nt)
    nelem = size(conn, 1)
    ngl = length(ω);  m = ngl - 2
    (m >= 1 && m^3 == no && size(K̂) == (ngl, ngl)) || return nothing
    Ki = K̂[2:ngl-1, 2:ngl-1];  wi = ω[2:ngl-1]
    h  = 1 ./ sqrt.(wi)
    E  = eigen(Symmetric(h .* Ki .* h'))
    S  = h .* E.vectors
    c  = zeros(3, nelem)
    ok = fill(true, nt)
    work = [_SCWork(size(A, 1), no, nb) for _ = 1:nt]
    lin(i, j, k) = i + (j - 1) * m + (k - 1) * m * m
    _sc_chunks(nelem, nt) do t, lo, hi
        w = work[t]
        for e = lo:hi
            ok[t] || return nothing
            _sc_extract_oo!(w, A, conn, e, nb, no, sgn)
            Aoo = w.Aoo
            if m == 1
                c[1, e] = Aoo[1, 1] / (Ki[1, 1] * wi[1]^2)
                Aoo[1, 1] > 0 || (ok[t] = false)
                continue
            end
            cx = Aoo[lin(2,1,1), lin(1,1,1)] / (Ki[2, 1] * wi[1] * wi[1])
            cy = Aoo[lin(1,2,1), lin(1,1,1)] / (Ki[2, 1] * wi[1] * wi[1])
            cz = Aoo[lin(1,1,2), lin(1,1,1)] / (Ki[2, 1] * wi[1] * wi[1])
            c[1, e] = cx;  c[2, e] = cy;  c[3, e] = cz
            (cx > 0 && cy > 0 && cz > 0) || (ok[t] = false;  continue)
            tol = 1e-10 * maximum(abs, Aoo)
            @inbounds for k2 = 1:m, j2 = 1:m, i2 = 1:m, k1 = 1:m, j1 = 1:m, i1 = 1:m
                v = 0.0
                j1 == j2 && k1 == k2 && (v += cx * Ki[i1, i2] * wi[j1] * wi[k1])
                i1 == i2 && k1 == k2 && (v += cy * wi[i1] * Ki[j1, j2] * wi[k1])
                i1 == i2 && j1 == j2 && (v += cz * wi[i1] * wi[j1] * Ki[k1, k2])
                if abs(Aoo[lin(i1, j1, k1), lin(i2, j2, k2)] - v) > tol
                    ok[t] = false;  break
                end
            end
        end
        return nothing
    end
    all(ok) || return nothing
    return _SCFDM(m, S, Matrix(S'), E.values, c)
end

# sgn·A_oo of element e -> w.Aoo
@inline function _sc_extract_oo!(w::_SCWork, A, conn, e, nb, no, sgn)
    @inbounds for j = 1:no;  w.loc[conn[e, nb+j]] = j;  end
    fill!(w.Aoo, 0.0)
    @inbounds for lc = 1:no, p in nzrange(A, conn[e, nb+lc])
        lr = w.loc[A.rowval[p]]
        lr > 0 && (w.Aoo[lr, lc] = sgn * A.nzval[p])
    end
    @inbounds for j = 1:no;  w.loc[conn[e, nb+j]] = 0;  end
    return nothing
end

# y[cols] = sgn·(A x)[cols] for a symmetric A (column dot products: threaded,
# no write conflicts); x must vanish on the fixed nodes
function _sc_colmatvec!(y, A::SparseMatrixCSC, x, cols::Vector{Int}, sgn, nt)
    _sc_chunks(length(cols), nt) do _, lo, hi
        @inbounds for q = lo:hi
            j = cols[q];  s = 0.0
            for p in nzrange(A, j);  s += A.nzval[p] * x[A.rowval[p]];  end
            y[j] = sgn * s
        end
        return nothing
    end
    return y
end

# per element: X[interior] ← α · A_oo⁻¹ Y[interior]   (threaded; interiors are disjoint)
function _sc_fdm_all!(X, Y, α, F::_SCFDM, fw::Vector{_SCFDMWork}, conn, nb, no, nt)
    nelem = size(conn, 1)
    _sc_chunks(nelem, min(nt, nelem)) do c, lo, hi
        w = fw[c]
        @inbounds for e = lo:hi
            for i = 1:no;  w.x[i] = Y[conn[e, nb+i]];  end
            _sc_fdm_solve!(w.x, w.t, F, e)
            for i = 1:no;  X[conn[e, nb+i]] = α * w.x[i];  end
        end
        return nothing
    end
    return X
end

"""
    el_sc_schur_cg!(u, A, f, conn, nb, ∂O, Γ, gΓ, K̂, ω; amg_method = "sa",
                    amg_rtol = 1e-12, itmax = 1000) -> Bool

Static condensation with AMG-CG, leveraging the SEM tensor product (see the
section header): CG on the skeleton Schur complement, applied matrix-free
with fast-diagonalization interior solves, preconditioned by the full-system
AMG V-cycle restricted to the skeleton. `K̂`, `ω`: 1D GLL stiffness and
weights. Returns false, doing nothing, unless every element is an affine
hexahedron of that tensor form (then use el_sc_solve!). Records :sc_extract
(checks, c's), :sc_condense (condensed right-hand side), :sc_factor (AMG
hierarchy), :sc_skeleton (CG) and :sc_recover; the CG iterations in
JX_AMG_STATS.
"""
function el_sc_schur_cg!(u::AbstractArray, A::SparseMatrixCSC{Float64}, f::AbstractVector,
                         conn::AbstractMatrix{<:Integer}, nb::Int,
                         ∂O::AbstractVector{<:Integer}, Γ::AbstractVector{<:Integer},
                         gΓ::AbstractVector, K̂::AbstractMatrix, ω::AbstractVector;
                         amg_method = "sa", amg_rtol::Real = 1e-12, itmax::Int = 1000)
    t0 = time_ns()
    nt = Threads.nthreads()
    R  = _sc_roles(A, conn, nb, ∂O, Γ, gΓ, nt)
    R === nothing && return false
    (; n, nelem, no, unk, ns, role, ug, sgn) = R
    nblas = BLAS.get_num_threads()
    F = _sc_fdm_setup(A, conn, nb, no, sgn, Matrix{Float64}(K̂), Vector{Float64}(ω), nt)
    F === nothing && return false
    intr = Int[conn[e, nb+i] for e = 1:nelem for i = 1:no]
    fw   = [_SCFDMWork(zeros(no), zeros(no), zeros(0, 0), Int32[]) for _ = 1:nt]
    JX_TIMINGS[:sc_extract] = (time_ns() - t0) / 1e9

    # ── condensed right-hand side ─────────────────────────────────────────────
    #   t = sgn (f − A_·Γ g);   rhs = t_b − Â_bo Â_oo⁻¹ t_o
    t0 = time_ns()
    t  = zeros(n);  xt = zeros(n);  y = zeros(n)
    @inbounds for i = 1:n;  role[i] == 2 || (t[i] = sgn * f[i]);  end
    @inbounds for (i, γ) in enumerate(Γ)
        gv = gΓ[i];  iszero(gv) && continue
        for p in nzrange(A, γ)
            r = A.rowval[p];  role[r] == 2 || (t[r] -= sgn * A.nzval[p] * gv)
        end
    end
    _sc_fdm_all!(xt, t, 1.0, F, fw, conn, nb, no, nt)          # xt_o = Â_oo⁻¹ t_o, 0 elsewhere
    _sc_colmatvec!(y, A, xt, unk, sgn, nt)                     # y_b = Â_bo xt_o
    rhs = Float64[t[unk[k]] - y[unk[k]] for k = 1:ns]
    JX_TIMINGS[:sc_condense] = (time_ns() - t0) / 1e9

    # ── AMG of the full system without the fixed nodes ───────────────────────
    keep = findall(!=(Int8(2)), role)
    kpos = zeros(Int, n);  kpos[keep] = 1:length(keep)
    P = jx_phase(:sc_factor) do
        Ak = A[keep, keep]
        sgn < 0 && (Ak.nzval .*= -1)
        jx_amg_setup(Ak; method = amg_method)
    end
    ub = jx_phase(:sc_skeleton) do
        nt > 1 && BLAS.set_num_threads(1)
        try
            _sc_pcg(ns, rhs, unk, kpos, length(keep), A, sgn, P, F, fw, conn, nb, no, intr,
                    xt, y, Float64(amg_rtol), itmax, nt)
        finally
            BLAS.set_num_threads(nblas)
        end
    end

    # ── interior recovery: u_o = Â_oo⁻¹ (t_o − Â_ob u_b) ───────────────────────
    t0 = time_ns()
    fill!(xt, 0.0)
    @inbounds for k = 1:ns;  xt[unk[k]] = ub[k];  ug[unk[k]] = ub[k];  end
    _sc_colmatvec!(y, A, xt, intr, sgn, nt)                    # y_o = Â_ob u_b
    @inbounds for g in intr;  y[g] = t[g] - y[g];  end
    _sc_fdm_all!(ug, y, 1.0, F, fw, conn, nb, no, nt)
    @inbounds for i = 1:n;  u[i] = ug[i];  end
    JX_TIMINGS[:sc_recover] = (time_ns() - t0) / 1e9
    return true
end

# Preconditioned CG on the skeleton (the stopping test of Krylov.cg:
# sqrt(rᵀ M⁻¹ r) ≤ rtol · its initial value)
function _sc_pcg(ns, b, unk, kpos, nk, A, sgn, P, F, fw, conn, nb, no, intr, xt, y, rtol, itmax, nt)
    x = zeros(ns);  r = copy(b);  z = zeros(ns);  p = zeros(ns);  q = zeros(ns)
    rk = zeros(nk);  zk = zeros(nk)
    precond! = (z, r) -> begin
        fill!(rk, 0.0)
        @inbounds for k = 1:ns;  rk[kpos[unk[k]]] = r[k];  end
        LinearAlgebra.ldiv!(zk, P.P, rk)
        @inbounds for k = 1:ns;  z[k] = zk[kpos[unk[k]]];  end
        z
    end
    schur! = (q, p) -> begin                       # q = S p
        fill!(xt, 0.0)
        @inbounds for k = 1:ns;  xt[unk[k]] = p[k];  end
        _sc_colmatvec!(y, A, xt, intr, sgn, nt)               # y_o = Â_ob p
        _sc_fdm_all!(xt, y, -1.0, F, fw, conn, nb, no, nt)    # xt_o = −Â_oo⁻¹ y_o
        _sc_colmatvec!(y, A, xt, unk, sgn, nt)                # y_b = Â_bb p + Â_bo xt_o
        @inbounds for k = 1:ns;  q[k] = y[unk[k]];  end
        q
    end
    precond!(z, r);  p .= z
    γ = dot(r, z);  γ0 = γ;  it = 0
    while sqrt(γ) > rtol * sqrt(γ0) && it < itmax
        schur!(q, p)
        α = γ / dot(p, q)
        axpy!(α, p, x);  axpy!(-α, q, r)
        precond!(z, r)
        γn = dot(r, z)
        β = γn / γ;  γ = γn
        @inbounds for k = 1:ns;  p[k] = z[k] + β * p[k];  end
        it += 1
    end
    rel = γ0 > 0 ? sqrt(γ / γ0) : 0.0
    JX_AMG_STATS[] = (iters = it, rel_resid = rel, levels = length(P.ml.levels) + 1,
                      method = Symbol(string(P.method, "+schur")))
    rel <= rtol || @warn " # Schur-complement CG did not converge to rtol=$rtol in $itmax iterations " *
                         "(final relative residual $rel)."
    return x
end
