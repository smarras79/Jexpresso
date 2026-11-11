#----------------------------------------------------------------------------------------------
# Dirichlet constraints detected from user_bc_dirichlet!  (1D, 2D, 3D; CPU, CG, non-Laguerre)
#
# No user input is needed: what a boundary tag constrains is inferred from what the case's
# user_bc_dirichlet! does to a node value.
#
#   * A component is CONSTRAINED on a boundary node if user_bc_dirichlet! changes its value.
#     The test is done once at setup on probe states (the initial state plus distinct offsets),
#     so a constraint that the current state already satisfies (e.g. w = 0 at rest) is still
#     detected. Values are still taken from user_bc_dirichlet! at every stage, so time-dependent
#     values keep working.
#   * NO-PENETRATION (free slip, 2D and 3D) is recognised when the returned (total) momentum
#     equals the tangential projection of the input momentum, m - (m.n) n, with n the boundary
#     normal. This covers axis-aligned walls and walls with topography. Such nodes get a node
#     normal built from all incident constrained boundary entities (sum of n n^T, eigen-
#     decomposed: nearly parallel normals merge into one averaged normal, distinct walls at
#     corners/edges stay separate), and the kernel itself projects the momentum and its
#     tendency with it. Tangential momentum is never frozen. In 1D a wall is ρu = 0, which the
#     component rule above already handles.
#   * Tags whose name contains "robin" (any case) are value+flux conditions: their values are
#     applied but no tendency is constrained, so a Neumann flux on the same tag acts.
#
# Tendency constraints are applied once per RHS evaluation AFTER the global assembly
# (DSS_global_RHS!), on masks that are themselves assembled across MPI ranks and periodic
# copies, so every copy of a node is treated identically.
#
# Momentum is assumed to be equations 2:(1+d) (CompEuler / shallow-water convention, d = 2 or
# 3); the no-penetration test is only attempted when user_bc_dirichlet! sets all of them, so
# other equation sets fall back to the component-wise rule.
#
# The boundary nodes and their tags/normals are collected once at setup, so the per-stage code
# is the same in every dimension; only the call to user_bc_dirichlet! (whose signature differs
# between 1D, 2D and 3D) is dispatched on the dimension.
#----------------------------------------------------------------------------------------------

const _BC_SENTINEL   = 4325789.0
const _BC_NP_SPLIT   = 0.2        # eigenvalue ratio separating distinct wall normals (~48 deg)

const _BC_PERIODIC_TAGS = ("periodicx", "periodic1", "periodicy", "periodic3", "periodicz", "periodic2")

bc_is_periodic_tag(tag) = tag in _BC_PERIODIC_TAGS
bc_is_robin_tag(tag)    = occursin("robin", lowercase(String(tag)))
bc_is_skipped_tag(tag)  = bc_is_periodic_tag(tag) || tag == "Laguerre"

struct St_bc_constraints{TF, SDT}
    SD      :: SDT                # NSD_1D(), NSD_2D() or NSD_3D()
    nd      :: Int                # momentum components handled by the projector (0, 2 or 3)
    # boundary node entries (one per (boundary entity, node); a node on several entities repeats)
    b_ip    :: Vector{Int}
    b_tag   :: Vector{String}
    b_n     :: Matrix{TF}         # (3, nb) outward normal; unused components are 0
    b_robin :: BitVector
    fixmask :: BitMatrix          # (npoin, neqs): value constrained -> tendency zeroed
    fix_ip  :: Vector{Int}        # flat list of the constrained (ip, ieq) pairs
    fix_ieq :: Vector{Int}
    np_ip   :: Vector{Int}        # no-penetration nodes
    np_P    :: Matrix{TF}         # (nd*nd, nnp) tangential projector I - sum v v^T, column-major
    is_np   :: BitVector          # (npoin)
end

#
# The user's routine has a different signature in each dimension.
#
@inline _bc_user!(::NSD_3D, qin, c, t, tag, qbdy, n1, n2, n3, xmin, xmax, ymin, ymax, zmin, zmax, qe_ip, SOL) =
    user_bc_dirichlet!(qin, c, t, tag, qbdy, n1, n2, n3, xmin, xmax, ymin, ymax, zmin, zmax, qe_ip, SOL)
@inline _bc_user!(::NSD_2D, qin, c, t, tag, qbdy, n1, n2, n3, xmin, xmax, ymin, ymax, zmin, zmax, qe_ip, SOL) =
    user_bc_dirichlet!(qin, c, t, tag, qbdy, n1, n2, qe_ip, SOL)
@inline _bc_user!(::NSD_1D, qin, c, t, tag, qbdy, n1, n2, n3, xmin, xmax, ymin, ymax, zmin, zmax, qe_ip, SOL) =
    user_bc_dirichlet!(qin, c, t, tag, qbdy, qe_ip, SOL)

@inline _bc_isset(v) = !(v == _BC_SENTINEL)

#
# Boundary node entries (ip, tag, normal) for each dimension.
#
function _bc_boundary_nodes(::NSD_3D, mesh, metrics, TF)
    ips = Int[]; tags = String[]; ns = NTuple{3,TF}[]
    for iface = 1:size(mesh.bdy_face_type, 1)
        tag = mesh.bdy_face_type[iface]
        bc_is_skipped_tag(tag) && continue
        for i = 1:mesh.ngl, j = 1:mesh.ngl
            push!(ips, mesh.poin_in_bdy_face[iface, i, j]); push!(tags, tag)
            push!(ns, (metrics.nx[iface, i, j], metrics.ny[iface, i, j], metrics.nz[iface, i, j]))
        end
    end
    return ips, tags, ns
end

function _bc_boundary_nodes(::NSD_2D, mesh, metrics, TF)
    ips = Int[]; tags = String[]; ns = NTuple{3,TF}[]
    for iedge = 1:mesh.nedges_bdy
        tag = mesh.bdy_edge_type[iedge]
        bc_is_skipped_tag(tag) && continue
        for k = 1:mesh.ngl
            push!(ips, mesh.poin_in_bdy_edge[iedge, k]); push!(tags, tag)
            push!(ns, (metrics.nx[iedge, k], metrics.ny[iedge, k], zero(TF)))
        end
    end
    return ips, tags, ns
end

function _bc_boundary_nodes(::NSD_1D, mesh, metrics, TF)
    return [1, mesh.npoin_linear], ["left", "right"], [(-one(TF), zero(TF), zero(TF)), (one(TF), zero(TF), zero(TF))]
end

_bc_np_dim(::NSD_3D) = 3
_bc_np_dim(::NSD_2D) = 2
_bc_np_dim(::NSD_1D) = 0

function build_bc_constraints(mesh, metrics, qn, qe, inputs, neqs,
                              xmin, xmax, ymin, ymax, zmin, zmax; t = 0.0)

    TF    = eltype(qe)
    SD    = mesh.SD
    npoin = mesh.npoin
    SOL   = inputs[:SOL_VARS_TYPE]
    lpert = (SOL == PERT())
    d     = _bc_np_dim(SD)
    lmom  = d > 0 && neqs >= 1 + d
    nS    = d*(d+1) ÷ 2               # independent entries of the symmetric S = sum n n^T

    b_ip, b_tag, b_nt = _bc_boundary_nodes(SD, mesh, metrics, TF)
    nb      = length(b_ip)
    b_n     = zeros(TF, 3, nb)
    b_robin = falses(nb)
    for e = 1:nb
        b_n[1, e], b_n[2, e], b_n[3, e] = b_nt[e]
        b_robin[e] = bc_is_robin_tag(b_tag[e])
    end

    ncol = neqs + nS + 1              # [component flags (neqs) | S entries | np count]
    W    = zeros(TF, npoin, ncol)

    nq    = size(qn, 2)
    qbdy  = zeros(TF, max(neqs, 7))
    probe = zeros(TF, nq)
    inval = zeros(TF, nq)
    m = zeros(TF, 3); o = zeros(TF, 3); b = zeros(TF, 3)
    # offset patterns: distinct magnitudes and signs so no constraint is met by accident
    pat = ((1.0, 1.37, 1.91, 2.53, 3.11, 3.79, 4.43, 5.17),
           (1.0, -1.61, 1.23, -2.07, 2.71, -3.29, 3.97, -4.61))

    for e = 1:nb
        b_robin[e] && continue
        ip  = b_ip[e]; tag = b_tag[e]
        n1, n2, n3 = b_n[1, e], b_n[2, e], b_n[3, e]
        nv  = (n1, n2, n3)
        qe_ip = @view(qe[ip, :])
        c_ip  = @view(mesh.coords[:, ip])

        changed  = falses(neqs)
        np_ok    = lmom
        np_moved = false

        for p = 0:length(pat)                 # p = 0: actual state, p >= 1: probes
            for k = 1:nq
                probe[k] = qn[ip, k]
            end
            if p > 0
                for k = 1:neqs
                    probe[k] += 1e-3 * pat[p][mod1(k, 8)] * (1 + abs(qn[ip, k]) + abs(qe[ip, k]))
                end
            end
            inval .= probe
            fill!(qbdy, _BC_SENTINEL)
            _bc_user!(SD, probe, c_ip, t, tag, qbdy, n1, n2, n3, xmin, xmax, ymin, ymax, zmin, zmax, qe_ip, SOL)

            for k = 1:neqs
                _bc_isset(qbdy[k]) || continue
                ref = abs(inval[k]) + abs(qe[ip, k]) + 1e-300
                abs(qbdy[k] - inval[k]) > 1e-10 * ref && (changed[k] = true)
            end

            # no-penetration test on the total momentum, components 2:(1+d)
            if np_ok
                if !all(c -> _bc_isset(qbdy[1+c]), 1:d)
                    np_ok = false
                elseif p > 0
                    mn = zero(TF); scale = zero(TF)
                    for c = 1:d
                        b[c] = lpert ? qe[ip, 1+c] : zero(TF)
                        m[c] = inval[1+c] + b[c]
                        o[c] = qbdy[1+c] + b[c]
                        mn  += m[c]*nv[c]; scale += abs(m[c])
                    end
                    scale += 1e-300
                    err = zero(TF)
                    for c = 1:d
                        err += abs(o[c] - (m[c] - mn*nv[c]))
                    end
                    if err > 1e-10 * scale
                        np_ok = false
                    elseif abs(mn) > 1e-10 * scale
                        np_moved = true
                    end
                end
            end
        end

        if np_ok && np_moved
            col = neqs
            for r = 1:d, c = r:d
                col += 1
                W[ip, col] += nv[r]*nv[c]
            end
            W[ip, neqs+nS+1] += 1
            for k = 1:neqs
                (2 <= k <= 1+d) && continue        # momentum handled by the projector
                changed[k] && (W[ip, k] += 1)
            end
        else
            for k = 1:neqs
                changed[k] && (W[ip, k] += 1)
            end
        end
    end

    # make every copy of a node (MPI ghost, periodic image) agree
    if SD != NSD_1D()
        cache = setup_assembler(SD, W, mesh.ip2gip, mesh.gip2owner)
        cache === nothing || assemble_mpi!(W, cache)
    end

    fixmask = falses(npoin, neqs)
    fix_ip  = Int[]; fix_ieq = Int[]
    for ip = 1:npoin, k = 1:neqs
        if W[ip, k] > 0.5
            fixmask[ip, k] = true
            push!(fix_ip, ip); push!(fix_ieq, k)
        end
    end

    is_np = falses(npoin)
    np_ip = Int[]
    Pcols = Vector{Vector{TF}}()
    if d > 0
        S = zeros(TF, d, d)
        for ip = 1:npoin
            W[ip, neqs+nS+1] > 0.5 || continue
            col = neqs
            for r = 1:d, c = r:d
                col += 1
                S[r, c] = S[c, r] = W[ip, col]
            end
            E = eigen(Symmetric(S))
            λmax = maximum(E.values)
            P = Matrix{TF}(I, d, d)
            for kk = 1:d
                if E.values[kk] > _BC_NP_SPLIT * λmax
                    v = E.vectors[:, kk]
                    # remove eigen-solver roundoff so flat walls get an exact axis normal
                    for c = 1:d
                        abs(v[c]) < 1e-12 && (v[c] = zero(TF))
                    end
                    v ./= sqrt(sum(abs2, v))
                    P .-= v * v'
                end
            end
            is_np[ip] = true
            push!(np_ip, ip)
            push!(Pcols, vec(P))
        end
    end
    np_P = zeros(TF, d*d, length(np_ip))
    for k = 1:length(np_ip)
        np_P[:, k] .= Pcols[k]
    end

    return St_bc_constraints{TF, typeof(SD)}(SD, d, b_ip, b_tag, b_n, b_robin,
                                              fixmask, fix_ip, fix_ieq, np_ip, np_P, is_np)
end

#
# Per RHS evaluation, before the RHS is computed: impose the constrained values on the state.
#
function apply_bc_constraints_state!(bc::St_bc_constraints, u, uaux, qe, t, coords,
                                     xmin, xmax, ymin, ymax, zmin, zmax, qbdy, neqs, npoin, inputs)
    SOL   = inputs[:SOL_VARS_TYPE]
    lpert = (SOL == PERT())
    SD    = bc.SD
    d     = bc.nd

    @inbounds for e = 1:length(bc.b_ip)
        ip = bc.b_ip[e]
        fill!(qbdy, _BC_SENTINEL)
        _bc_user!(SD, @view(uaux[ip, :]), @view(coords[:, ip]), t, bc.b_tag[e], qbdy,
                  bc.b_n[1, e], bc.b_n[2, e], bc.b_n[3, e],
                  xmin, xmax, ymin, ymax, zmin, zmax, @view(qe[ip, :]), SOL)
        robin = bc.b_robin[e]
        for k = 1:neqs
            _bc_isset(qbdy[k]) || continue
            if robin || bc.fixmask[ip, k]
                uaux[ip, k] = qbdy[k]
            end
        end
    end

    # no-penetration: project the total momentum with the node projector
    if d == 3
        @inbounds for kk = 1:length(bc.np_ip)
            ip = bc.np_ip[kk]
            P  = @view bc.np_P[:, kk]
            b1 = lpert ? qe[ip, 2] : zero(eltype(uaux))
            b2 = lpert ? qe[ip, 3] : zero(eltype(uaux))
            b3 = lpert ? qe[ip, 4] : zero(eltype(uaux))
            m1 = uaux[ip, 2] + b1; m2 = uaux[ip, 3] + b2; m3 = uaux[ip, 4] + b3
            uaux[ip, 2] = P[1]*m1 + P[4]*m2 + P[7]*m3 - b1
            uaux[ip, 3] = P[2]*m1 + P[5]*m2 + P[8]*m3 - b2
            uaux[ip, 4] = P[3]*m1 + P[6]*m2 + P[9]*m3 - b3
        end
    elseif d == 2
        @inbounds for kk = 1:length(bc.np_ip)
            ip = bc.np_ip[kk]
            P  = @view bc.np_P[:, kk]
            b1 = lpert ? qe[ip, 2] : zero(eltype(uaux))
            b2 = lpert ? qe[ip, 3] : zero(eltype(uaux))
            m1 = uaux[ip, 2] + b1; m2 = uaux[ip, 3] + b2
            uaux[ip, 2] = P[1]*m1 + P[3]*m2 - b1
            uaux[ip, 3] = P[2]*m1 + P[4]*m2 - b2
        end
    end

    uaux2u!(u, uaux, neqs, npoin)
    return nothing
end

#
# Per RHS evaluation, after DSS_global_RHS! (and before the mass division): constrain the tendency.
#
function apply_bc_constraints_rhs!(RHS, bc::St_bc_constraints)
    d = bc.nd
    if d == 3
        @inbounds for kk = 1:length(bc.np_ip)
            ip = bc.np_ip[kk]
            P  = @view bc.np_P[:, kk]
            r1 = RHS[ip, 2]; r2 = RHS[ip, 3]; r3 = RHS[ip, 4]
            RHS[ip, 2] = P[1]*r1 + P[4]*r2 + P[7]*r3
            RHS[ip, 3] = P[2]*r1 + P[5]*r2 + P[8]*r3
            RHS[ip, 4] = P[3]*r1 + P[6]*r2 + P[9]*r3
        end
    elseif d == 2
        @inbounds for kk = 1:length(bc.np_ip)
            ip = bc.np_ip[kk]
            P  = @view bc.np_P[:, kk]
            r1 = RHS[ip, 2]; r2 = RHS[ip, 3]
            RHS[ip, 2] = P[1]*r1 + P[3]*r2
            RHS[ip, 3] = P[2]*r1 + P[4]*r2
        end
    end
    @inbounds for kk = 1:length(bc.fix_ip)
        RHS[bc.fix_ip[kk], bc.fix_ieq[kk]] = zero(eltype(RHS))
    end
    return nothing
end
