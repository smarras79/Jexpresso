#----------------------------------------------------------------------------------------------
# Dirichlet constraints detected from user_bc_dirichlet!  (3D, CPU, non-Laguerre)
#
# No user input is needed: what a boundary tag constrains is inferred from what the case's
# user_bc_dirichlet! does to a node value.
#
#   * A component is CONSTRAINED on a boundary node if user_bc_dirichlet! changes its value.
#     The test is done once at setup on probe states (the initial state plus distinct offsets),
#     so a constraint that the current state already satisfies (e.g. w = 0 at rest) is still
#     detected. Values are still taken from user_bc_dirichlet! at every stage, so time-dependent
#     values keep working.
#   * NO-PENETRATION (free slip) is recognised when the returned (total) momentum equals the
#     tangential projection of the input momentum, m - (m.n) n, with n the face normal. This
#     covers axis-aligned walls and walls with topography. Such nodes get a node normal built
#     from all incident constrained faces (sum of n n^T, eigen-decomposed: nearly parallel face
#     normals merge into one averaged normal, distinct walls at edges/corners stay separate),
#     and the kernel itself projects the momentum and its tendency with it. Tangential
#     momentum is never frozen.
#   * Tags whose name contains "robin" (any case) are value+flux conditions: their values are
#     applied but no tendency is constrained, so a Neumann flux on the same tag acts.
#
# Tendency constraints are applied once per RHS evaluation AFTER the global assembly
# (DSS_global_RHS!), on masks that are themselves assembled across MPI ranks and periodic
# copies, so every copy of a node is treated identically.
#
# Momentum is assumed to be equations 2:4 (CompEuler convention); the no-penetration test is only
# attempted when user_bc_dirichlet! sets all three of them, so other equation sets fall back to
# the component-wise rule.
#----------------------------------------------------------------------------------------------

const _BC_SENTINEL   = 4325789.0
const _BC_NP_SPLIT   = 0.2        # eigenvalue ratio separating distinct wall normals (~48 deg)

const _BC_PERIODIC_TAGS = ("periodicx", "periodic1", "periodicy", "periodic3", "periodicz", "periodic2")

bc_is_periodic_tag(tag) = tag in _BC_PERIODIC_TAGS
# case-insensitive "robin" substring test without allocating (called per boundary face and RK stage)
function bc_is_robin_tag(tag)
    s = String(tag)                       # no copy when tag is already a String
    n = ncodeunits(s)
    @inbounds for i in 1:n-4
        match = true
        for k in 1:5
            if lowercase(Char(codeunit(s, i+k-1))) != ('r', 'o', 'b', 'i', 'n')[k]
                match = false
                break
            end
        end
        match && return true
    end
    return false
end

struct St_bc_constraints{TF}
    fixmask :: BitMatrix          # (npoin, neqs): value constrained -> tendency zeroed
    fix_ip  :: Vector{Int}        # flat list of the constrained (ip, ieq) pairs
    fix_ieq :: Vector{Int}
    np_ip   :: Vector{Int}        # no-penetration nodes
    np_P    :: Matrix{TF}         # (9, nnp) tangential projector I - sum v v^T, column-major 3x3
    is_np   :: BitVector          # (npoin)
end

#
# Probe one boundary node: call user_bc_dirichlet! on `qin`, return the qbdy it produced.
#
@inline function _bc_call!(qbdy, qin, coords_ip, t, tag, n1, n2, n3,
                           xmin, xmax, ymin, ymax, zmin, zmax, qe_ip, SOL_VARS_TYPE)
    fill!(qbdy, _BC_SENTINEL)
    user_bc_dirichlet!(qin, coords_ip, t, tag, qbdy, n1, n2, n3,
                       xmin, xmax, ymin, ymax, zmin, zmax, qe_ip, SOL_VARS_TYPE)
    return qbdy
end

@inline _bc_isset(v) = !(v == _BC_SENTINEL)

function build_bc_constraints(mesh, metrics, qn, qe, inputs, neqs,
                              xmin, xmax, ymin, ymax, zmin, zmax; t = 0.0)

    TF    = eltype(qe)
    npoin = mesh.npoin
    ngl   = mesh.ngl
    SOL   = inputs[:SOL_VARS_TYPE]
    lpert = (SOL == PERT())
    lmom  = neqs >= 4

    ncol = neqs + 7                       # [component flags (neqs) | S = sum n n^T (6) | np count]
    W    = zeros(TF, npoin, ncol)

    nq    = size(qn, 2)
    qbdy  = zeros(TF, max(neqs, 7))
    probe = zeros(TF, nq)
    # offset patterns: distinct magnitudes and signs so no constraint is met by accident
    pat = ((1.0, 1.37, 1.91, 2.53, 3.11, 3.79, 4.43, 5.17),
           (1.0, -1.61, 1.23, -2.07, 2.71, -3.29, 3.97, -4.61))

    for iface = 1:size(mesh.bdy_face_type, 1)
        tag = mesh.bdy_face_type[iface]
        (bc_is_periodic_tag(tag) || bc_is_robin_tag(tag)) && continue
        for i = 1:ngl, j = 1:ngl
            ip = mesh.poin_in_bdy_face[iface, i, j]
            n1 = metrics.nx[iface, i, j]; n2 = metrics.ny[iface, i, j]; n3 = metrics.nz[iface, i, j]
            qe_ip = @view(qe[ip, :])
            c_ip  = @view(mesh.coords[ip, :])

            changed  = falses(neqs)
            np_ok    = lmom
            np_moved = false

            for p = 0:length(pat)                 # p = 0: actual state, p >= 1: probes
                for k = 1:nq
                    probe[k] = qn[ip, k]
                end
                if p > 0
                    for k = 1:neqs
                        a = 1e-3 * pat[p][mod1(k, 8)] * (1 + abs(qn[ip, k]) + abs(qe[ip, k]))
                        probe[k] += a
                    end
                end
                inval = copy(probe)
                _bc_call!(qbdy, probe, c_ip, t, tag, n1, n2, n3,
                          xmin, xmax, ymin, ymax, zmin, zmax, qe_ip, SOL)

                for k = 1:neqs
                    _bc_isset(qbdy[k]) || continue
                    ref = abs(inval[k]) + abs(qe[ip, k]) + 1e-300
                    abs(qbdy[k] - inval[k]) > 1e-10 * ref && (changed[k] = true)
                end

                # no-penetration test on the total momentum
                if np_ok
                    if !(_bc_isset(qbdy[2]) && _bc_isset(qbdy[3]) && _bc_isset(qbdy[4]))
                        np_ok = false
                    elseif p > 0
                        b1 = lpert ? qe[ip, 2] : zero(TF)
                        b2 = lpert ? qe[ip, 3] : zero(TF)
                        b3 = lpert ? qe[ip, 4] : zero(TF)
                        m1 = inval[2] + b1; m2 = inval[3] + b2; m3 = inval[4] + b3
                        mn = m1*n1 + m2*n2 + m3*n3
                        e1 = m1 - mn*n1;  e2 = m2 - mn*n2;  e3 = m3 - mn*n3
                        o1 = qbdy[2] + b1; o2 = qbdy[3] + b2; o3 = qbdy[4] + b3
                        scale = abs(m1) + abs(m2) + abs(m3) + 1e-300
                        if abs(o1 - e1) + abs(o2 - e2) + abs(o3 - e3) > 1e-10 * scale
                            np_ok = false
                        elseif abs(mn) > 1e-10 * scale
                            np_moved = true
                        end
                    end
                end
            end

            if np_ok && np_moved
                W[ip, neqs+1] += n1*n1; W[ip, neqs+2] += n2*n2; W[ip, neqs+3] += n3*n3
                W[ip, neqs+4] += n1*n2; W[ip, neqs+5] += n1*n3; W[ip, neqs+6] += n2*n3
                W[ip, neqs+7] += 1
                for k = 1:neqs
                    (2 <= k <= 4) && continue           # momentum handled by the projector
                    changed[k] && (W[ip, k] += 1)
                end
            else
                for k = 1:neqs
                    changed[k] && (W[ip, k] += 1)
                end
            end
        end
    end

    # make every copy of a node (MPI ghost, periodic image) agree
    if mesh.SD != NSD_1D()
        cache = setup_assembler(mesh.SD, W, mesh.ip2gip, mesh.gip2owner)
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
    Pcols = Vector{NTuple{9,TF}}()
    S = zeros(TF, 3, 3)
    for ip = 1:npoin
        W[ip, neqs+7] > 0.5 || continue
        S[1,1] = W[ip, neqs+1]; S[2,2] = W[ip, neqs+2]; S[3,3] = W[ip, neqs+3]
        S[1,2] = S[2,1] = W[ip, neqs+4]
        S[1,3] = S[3,1] = W[ip, neqs+5]
        S[2,3] = S[3,2] = W[ip, neqs+6]
        E = eigen(Symmetric(S))
        λmax = maximum(E.values)
        P = Matrix{TF}(I, 3, 3)
        for kk = 1:3
            if E.values[kk] > _BC_NP_SPLIT * λmax
                v = E.vectors[:, kk]
                # remove eigen-solver roundoff so flat walls get an exact axis normal
                for c = 1:3
                    abs(v[c]) < 1e-12 && (v[c] = zero(TF))
                end
                v ./= sqrt(sum(abs2, v))
                P .-= v * v'
            end
        end
        is_np[ip] = true
        push!(np_ip, ip)
        push!(Pcols, Tuple(vec(P)))
    end
    np_P = zeros(TF, 9, length(np_ip))
    for k = 1:length(np_ip)
        np_P[:, k] .= Pcols[k]
    end

    return St_bc_constraints{TF}(fixmask, fix_ip, fix_ieq, np_ip, np_P, is_np)
end

#
# Per RHS evaluation, before the RHS is computed: impose the constrained values on the state.
#
function apply_bc_constraints_state!(bc::St_bc_constraints, u, uaux, qe, t, coords,
                                     nx, ny, nz, mesh_bdy_face_type, poin_in_bdy_face, ngl,
                                     xmin, xmax, ymin, ymax, zmin, zmax, qbdy, neqs, npoin, inputs)
    # function barrier: SOL_VARS_TYPE comes from the untyped inputs Dict
    _apply_bc_constraints_state!(bc, u, uaux, qe, t, coords,
                                 nx, ny, nz, mesh_bdy_face_type, poin_in_bdy_face, ngl,
                                 xmin, xmax, ymin, ymax, zmin, zmax, qbdy, neqs, npoin, inputs[:SOL_VARS_TYPE])
end

function _apply_bc_constraints_state!(bc::St_bc_constraints, u, uaux, qe, t, coords,
                                      nx, ny, nz, mesh_bdy_face_type, poin_in_bdy_face, ngl,
                                      xmin, xmax, ymin, ymax, zmin, zmax, qbdy, neqs, npoin, SOL)
    lpert = (SOL == PERT())

    for iface = 1:size(mesh_bdy_face_type, 1)
        tag = mesh_bdy_face_type[iface]
        bc_is_periodic_tag(tag) && continue
        robin = bc_is_robin_tag(tag)
        for i = 1:ngl, j = 1:ngl
            ip = poin_in_bdy_face[iface, i, j]
            fill!(qbdy, _BC_SENTINEL)
            user_bc_dirichlet!(@view(uaux[ip, :]), @view(coords[ip, :]), t, tag, qbdy,
                               nx[iface, i, j], ny[iface, i, j], nz[iface, i, j],
                               xmin, xmax, ymin, ymax, zmin, zmax, @view(qe[ip, :]), SOL)
            for k = 1:neqs
                _bc_isset(qbdy[k]) || continue
                if robin || bc.fixmask[ip, k]
                    uaux[ip, k] = qbdy[k]
                end
            end
        end
    end

    # no-penetration: project the total momentum with the node projector
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

    uaux2u!(u, uaux, neqs, npoin)
    return nothing
end

#
# Per RHS evaluation, after DSS_global_RHS! (and before the mass division): constrain the tendency.
#
function apply_bc_constraints_rhs!(RHS, bc::St_bc_constraints)
    @inbounds for kk = 1:length(bc.np_ip)
        ip = bc.np_ip[kk]
        P  = @view bc.np_P[:, kk]
        r1 = RHS[ip, 2]; r2 = RHS[ip, 3]; r3 = RHS[ip, 4]
        RHS[ip, 2] = P[1]*r1 + P[4]*r2 + P[7]*r3
        RHS[ip, 3] = P[2]*r1 + P[5]*r2 + P[8]*r3
        RHS[ip, 4] = P[3]*r1 + P[6]*r2 + P[9]*r3
    end
    @inbounds for kk = 1:length(bc.fix_ip)
        RHS[bc.fix_ip[kk], bc.fix_ieq[kk]] = zero(eltype(RHS))
    end
    return nothing
end
