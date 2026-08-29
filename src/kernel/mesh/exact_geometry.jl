#---------------------------------------------------------------------------------
# exact_geometry.jl
#
# Curve the high-order nodes of a 2D grid onto the EXACT geometry that the mesh
# generator only approximated by straight sides.
#
# WHY
# ---
# gmsh writes a linear (straight-sided) grid: a `Circle` arc in the .geo comes
# back as a polyline whose vertices happen to sit on the circle. Populating that
# element with LGL nodes puts every high-order node on the CHORD, so no matter
# how large `:nop` is, the boundary the solver sees is a polygon. The wall
# normal is then wrong by O(h) at every node and a slip wall spuriously
# generates vorticity at each polygon corner. This is the 2D analogue of what
# `project_nodes_to_shell!` already does for a spherical shell, where the
# high-order nodes are pushed back out onto the sphere.
#
# HOW  (Kopriva, J. Sci. Comput. 26(3):301-327, 2006)
# ---------------------------------------------------------------------------
# Section 3 of that paper states the standard element mapping: a LINEAR
# BLENDING TRANSFINITE map (Gordon & Hall, 1973) between the four edge curves,
# which is a polynomial in each direction when the curved sides are. That is
# what is built here, in two steps:
#
#   1. Put the N+1 nodes of a boundary EDGE onto the exact curve. The degree-N
#      Lagrange interpolant through them is the boundary representation Γ(ξ) —
#      isoparametric, the same order as the solution. Kopriva's Remark 2 is
#      explicit that this is enough and that representing the boundary to
#      HIGHER order than N does not help.
#
#   2. Fill the element interior by blending the four edge curves. Linear
#      blending of curves in P^N stays in P^N.
#
# Theorem 3 is what says this is safe: in 2D the cross-product metric terms
#
#       J a¹ = Y_η x̂ - X_η ŷ,      J a² = -Y_ξ x̂ + X_ξ ŷ                  (38)
#
# satisfy the discrete metric identities — and the scheme is therefore
# constant-state preserving — iff the mapping they are computed from is in P^N.
# `build_metric_terms!(…, NSD_2D)` differentiates the degree-N interpolant
# through `mesh.x/y`, so its mapping is in P^N however the nodes are placed,
# and curving the grid this way keeps free-stream preservation exactly. (What
# Remark 4 warns against is the other route: substituting metric terms obtained
# by analytically differentiating the true circular map. Jexpresso never does
# that, and this file must not start.)
#
# Written in terms of the DISPLACEMENT δ = (curved node) - (straight node),
# step 2 is
#
#   δ(ξ,η) = (1-η)/2 δ₁(ξ) + (1+η)/2 δ₂(ξ) + (1-ξ)/2 δ₃(η) + (1+ξ)/2 δ₄(η)
#            - [bilinear corner terms]
#
# and the corner terms VANISH here, because the element corners are the linear
# gmsh vertices and those are left exactly where gmsh put them (δ = 0 there).
# Two consequences, both of which are why this routine needs no communication
# and no neighbour bookkeeping:
#
#   * the blend of a curved edge is identically zero on the two edges adjacent
#     to it, so a neighbouring element that does NOT touch the curved boundary
#     sees its shared edge unmoved — the grid stays exactly conforming;
#   * an element with two curved edges is handled by the plain sum above, with
#     no double counting.
#
# Only the interior nodes of a boundary element, and the interior nodes of the
# curved boundary edges themselves, ever move. Everything else is untouched.
#
# INPUT
# -----
#   :exact_geometry => Dict(<boundary tag> => <shape spec>)
#
# where <boundary tag> is the gmsh `Physical Curve` name and <shape spec> is
# either
#
#   :circle                       centre and radius are FITTED (least squares)
#                                 from the linear vertices the mesh file
#                                 already carries on that boundary, and the fit
#                                 is rejected if they do not actually lie on a
#                                 circle — the same "trust the grid, but check
#                                 it" policy as project_nodes_to_shell!;
#
#   (:circle, xc, yc, r)          centre and radius given explicitly.
#
# Example (problems/CompEuler/shock_circle/user_inputs.jl):
#
#   :exact_geometry => Dict("circle_boundary" => :circle),
#
# No key, or an empty Dict, means "do nothing" — every existing case is
# unaffected.
#
# LIMITS
# ------
#   * 2D only. The 3D construction blends six FACE surfaces rather than four
#     edge curves, and Theorem 5 of the same paper shows the cross-product
#     metrics do NOT satisfy the metric identities there whatever the mapping
#     is — 3D needs the conservative curl form as well, so it is its own change.
#     NSD_1D and NSD_3D are explicit no-ops rather than silent ones.
#
#   * ONE ELEMENT LAYER. Only the elements that own a curved edge are curved,
#     which is the standard construction (it is the mesh of Kopriva's Fig. 2)
#     and is what keeps the correction local and the grid conforming. The grid
#     must therefore already resolve the curvature: an element thinner, normal
#     to the wall, than the sagitta of its own boundary arc folds, and
#     _check_curved_elements stops the run rather than hand a negative Jacobian
#     to build_metric_terms!.
#
#   * AMR. This runs inside mod_mesh_read_gmsh!, so an adapted grid is curved
#     again from its own boundary vertices each time it is rebuilt. Those
#     vertices come from refining the STRAIGHT-SIDED parent, so refinement
#     brings the vertices onto the wall at the usual O(h²) and the curving then
#     corrects the nodes between them — correct, but the h in "O(h^{N+1}) from
#     the true circle" is the parent element's, not the child's, until the
#     vertices themselves are put on the geometry. Untested with :ladapt.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Shapes.
#
# A shape only has to answer one question: given a point, where is the nearest
# point of the exact geometry? Adding :sphere / :cylinder for the 3D case is a
# matter of adding a struct and a `snap_to_shape` method.
#---------------------------------------------------------------------------------
struct ExactCircle{T<:AbstractFloat}
    xc::T
    yc::T
    r::T
end

Base.show(io::IO, c::ExactCircle) =
    print(io, "circle(centre = (", c.xc, ", ", c.yc, "), r = ", c.r, ")")

# Radial projection onto the circle. This is the 2D twin of the radial snap in
# project_nodes_to_shell!, and it is shape-generic: any geometry that can name
# its nearest point plugs in here unchanged. Theorem 3 does not care WHICH
# parametrization of the arc the N+1 nodes end up at — only that they lie on
# the curve, so that the interpolant through them is the P^N boundary.
function snap_to_shape(c::ExactCircle, x, y)
    dx = x - c.xc
    dy = y - c.yc
    d  = sqrt(dx*dx + dy*dy)
    d > 0 || return (x, y)              # node sits on the centre: nothing to do
    s = c.r/d
    return (c.xc + s*dx, c.yc + s*dy)
end


# Characteristic length of a shape, for turning an absolute distance into a
# relative one.
_shape_scale(c::ExactCircle) = c.r


#---------------------------------------------------------------------------------
# Least-squares circle through a point cloud (Kåsa / algebraic fit).
#
# Minimises Σ (x² + y² + Dx + Ey + F)², i.e. the 3x3 normal system
#
#   [ Σx²  Σxy  Σx ] [D]     [ Σ x(x²+y²) ]
#   [ Σxy  Σy²  Σy ] [E] = - [ Σ y(x²+y²) ]
#   [ Σx   Σy   n  ] [F]     [ Σ  (x²+y²) ]
#
# with xc = -D/2, yc = -E/2, r = sqrt(xc² + yc² - F). The sums are formed
# locally and Allreduced, so a partitioned grid fits the same circle on every
# rank and the snapped coordinates agree bit-for-bit across a partition
# boundary.
#
# Returns (circle, max relative radial residual) or `nothing` when the cloud is
# degenerate (fewer than 3 points, or collinear — a straight boundary, for
# which the normal system is singular).
#---------------------------------------------------------------------------------
function _fit_circle_mpi(xs::AbstractVector, ys::AbstractVector, comm)

    s = zeros(Float64, 7)               # n, Σx, Σy, Σx², Σxy, Σy², Σ(x²+y²)
    t = zeros(Float64, 3)               # Σx·q, Σy·q, Σq   with q = x²+y²
    for k in eachindex(xs)
        x = Float64(xs[k]); y = Float64(ys[k]); q = x*x + y*y
        s[1] += 1.0; s[2] += x;   s[3] += y
        s[4] += x*x; s[5] += x*y; s[6] += y*y; s[7] += q
        t[1] += x*q; t[2] += y*q; t[3] += q
    end
    s = MPI.Allreduce(s, MPI.SUM, comm)
    t = MPI.Allreduce(t, MPI.SUM, comm)

    n = s[1]
    n >= 3.0 || return nothing

    A = [s[4] s[5] s[2];
         s[5] s[6] s[3];
         s[2] s[3] n]
    b = -[t[1], t[2], t[3]]

    # Collinear (or nearly so) => singular A => not a circle. Scale the
    # determinant test by the cloud's own size so it is dimensionless.
    scale = max(s[4] + s[6] - (s[2]^2 + s[3]^2)/n, eps(Float64))   # n·variance
    abs(det(A)) > 1.0e-12 * n * scale^2 || return nothing

    D, E, F = A \ b
    xc = -0.5*D
    yc = -0.5*E
    disc = xc*xc + yc*yc - F
    disc > 0.0 || return nothing
    r = sqrt(disc)
    r > 0.0 || return nothing

    # Residual: how far the fitted circle actually is from the vertices it was
    # fitted to. A polygonal boundary that is not a circle shows up here.
    resid = 0.0
    for k in eachindex(xs)
        resid = max(resid, abs(sqrt((Float64(xs[k]) - xc)^2 + (Float64(ys[k]) - yc)^2) - r))
    end
    resid = MPI.Allreduce(resid, MPI.MAX, comm)

    return (ExactCircle{TFloat}(TFloat(xc), TFloat(yc), TFloat(r)), resid/r)
end


#---------------------------------------------------------------------------------
# D[k,i] = ℓ_i'(ξ_k): the nodal derivative matrix of the degree-N Lagrange basis
# on the interpolation nodes, in barycentric form (Berrut & Trefethen 2004,
# §9.3). Used here to check element validity, and by build_metric_terms! to get
# the tangent along a curved boundary edge.
#---------------------------------------------------------------------------------
function lagrange_nodal_derivative_matrix(ξ::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(ξ)
    w = ones(T, n)
    @inbounds for j = 1:n, k = 1:n
        k == j && continue
        w[j] /= (ξ[j] - ξ[k])
    end
    D = zeros(T, n, n)
    @inbounds for k = 1:n
        acc = zero(T)
        for i = 1:n
            i == k && continue
            D[k,i] = (w[i]/w[k])/(ξ[k] - ξ[i])
            acc += D[k,i]
        end
        D[k,k] = -acc                    # rows sum to zero: d/dξ of a constant
    end
    return D
end


#---------------------------------------------------------------------------------
# Resolve one <boundary tag> => <shape spec> entry into a concrete shape.
#---------------------------------------------------------------------------------
function _resolve_shape(tag::AbstractString, spec, xs, ys, comm, rank, suppress)

    # (:circle, xc, yc, r) — the case deck states the geometry outright.
    if spec isa Tuple && length(spec) == 4 && spec[1] === :circle
        return ExactCircle{TFloat}(TFloat(spec[2]), TFloat(spec[3]), TFloat(spec[4]))
    end

    # :circle — fit it from what the grid file itself says.
    if spec === :circle
        fit = _fit_circle_mpi(xs, ys, comm)
        if fit === nothing
            println_rank(string(" #   :exact_geometry \"", tag,
                                "\" => :circle ignored: its vertices are collinear or too few ",
                                "to define a circle.");
                         msg_rank = rank, suppress = suppress)
            return nothing
        end
        circle, rresid = fit
        if rresid > 1.0e-6
            println_rank(string(" #   :exact_geometry \"", tag,
                                "\" => :circle ignored: the boundary vertices are not on a circle ",
                                "(best fit ", circle, " leaves a relative residual of ", rresid, ").");
                         msg_rank = rank, suppress = suppress)
            return nothing
        end
        return circle
    end

    error(" # ERROR exact_geometry.jl: :exact_geometry[\"" * String(tag) * "\"] => " *
          string(spec) * " is not a shape I know.\n" *
          " #   Use :circle (centre and radius fitted from the grid) or\n" *
          " #   (:circle, xc, yc, r) to state them explicitly.")
end


#---------------------------------------------------------------------------------
# snap_nodes_to_exact_geometry!(mesh, lgl, inputs, SD)
#
# The entry point, called from mod_mesh_read_gmsh! once the high-order nodes,
# `connijk`, and the boundary-edge tables all exist. Rewrites `mesh.x/y` in
# place; `mesh.coords` is filled from x/y at the end of mod_mesh_mesh_driver,
# so the curved coordinates propagate everywhere (metrics, BCs, VTK) on their
# own.
#
# 1D and 3D are no-ops for now: the 3D construction blends six FACE surfaces
# rather than four edge curves, and Kopriva's Theorem 5 shows the cross-product
# metrics do NOT satisfy the metric identities there, so it needs the curl form
# as well and is deliberately left for its own change.
#---------------------------------------------------------------------------------
snap_nodes_to_exact_geometry!(mesh::St_mesh, lgl, inputs::Dict{Symbol,Any}, ::NSD_1D) = nothing
snap_nodes_to_exact_geometry!(mesh::St_mesh, lgl, inputs::Dict{Symbol,Any}, ::NSD_3D) = nothing

function snap_nodes_to_exact_geometry!(mesh::St_mesh, lgl, inputs::Dict{Symbol,Any}, ::NSD_2D)

    spec = get(inputs, :exact_geometry, nothing)
    (spec === nothing || isempty(spec)) && return nothing

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    ngl = mesh.ngl
    if ngl < 3
        # nop = 1: the grid has no nodes other than the gmsh vertices, which
        # already lie on the geometry. Nothing to curve.
        println_rank(" #   :exact_geometry ignored at :nop => 1 (no high-order nodes to move).";
                     msg_rank = rank, suppress = mesh.msg_suppress)
        return nothing
    end

    println_rank(" # SNAP HIGH-ORDER NODES ONTO EXACT GEOMETRY ....................";
                 msg_rank = rank, suppress = mesh.msg_suppress)

    #-----------------------------------------------------------------------------
    # Displacement field. Nonzero ONLY on the interior nodes of a curved
    # boundary edge; the gmsh vertices keep their coordinates so that δ = 0 at
    # every element corner, which is what makes the transfinite blend below
    # conforming (see the header).
    #-----------------------------------------------------------------------------
    δx = zeros(TFloat, mesh.npoin)
    δy = zeros(TFloat, mesh.npoin)
    on_curved = falses(mesh.npoin)

    # Sorted by tag, so every rank walks them in the same order: the circle fit
    # below is an MPI collective and must be entered in lockstep. Keys may be
    # written as strings or symbols in the deck; both name the same gmsh group.
    ncurved_tags = 0
    entries = sort!([(String(k), v) for (k, v) in spec]; by = first)
    for (tagstr, shape_spec) in entries

        edges = findall(t -> t !== nothing && String(t) == tagstr, mesh.bdy_edge_type)

        # Vertices of this boundary, for the fit. Only the LINEAR ones: they are
        # the geometry the mesh file asserts; the high-order ones are the chord
        # interpolant we are about to correct.
        xs = TFloat[]; ys = TFloat[]
        for iedge in edges
            for igl in (1, ngl)
                ip = mesh.poin_in_bdy_edge[iedge, igl]
                push!(xs, mesh.x[ip]); push!(ys, mesh.y[ip])
            end
        end

        shape = _resolve_shape(tagstr, shape_spec, xs, ys, comm, rank, mesh.msg_suppress)
        shape === nothing && continue

        # A tag that exists in no .msh anywhere is a typo in the deck, not a
        # partitioning artifact — check it globally before reporting.
        nedges_glob = MPI.Allreduce(length(edges), MPI.SUM, comm)
        if nedges_glob == 0
            println_rank(string(" #   :exact_geometry \"", tagstr,
                                "\" ignored: no boundary edge carries that tag.");
                         msg_rank = rank, suppress = mesh.msg_suppress)
            continue
        end
        ncurved_tags += 1

        dmax = 0.0
        dvert = 0.0
        for iedge in edges
            for igl = 1:ngl
                ip = mesh.poin_in_bdy_edge[iedge, igl]
                on_curved[ip] = true
                xs_new, ys_new = snap_to_shape(shape, mesh.x[ip], mesh.y[ip])
                if igl == 1 || igl == ngl
                    # THE LINEAR VERTICES DO NOT MOVE. They are the element
                    # corners, and the blend below is conforming precisely
                    # because δ vanishes there (see the header): snapping a
                    # corner would revive the Gordon-Hall corner terms and pull
                    # the shared side edges of neighbouring elements apart.
                    # gmsh already puts them on the geometry, so there is
                    # nothing to gain — but say so if this grid does not.
                    dvert = max(dvert, hypot(xs_new - mesh.x[ip], ys_new - mesh.y[ip]))
                    continue
                end
                δx[ip] = xs_new - mesh.x[ip]
                δy[ip] = ys_new - mesh.y[ip]
                dmax = max(dmax, hypot(δx[ip], δy[ip]))
            end
        end

        println_rank(string(" #   \"", tagstr, "\" -> ", shape,
                            " ; ", nedges_glob, " edges ; largest node correction = ",
                            MPI.Allreduce(dmax, MPI.MAX, comm));
                     msg_rank = rank, suppress = mesh.msg_suppress)

        dvert = MPI.Allreduce(dvert, MPI.MAX, comm)
        if dvert > 1.0e-8*_shape_scale(shape)
            println_rank(string(" #   NOTE \"", tagstr, "\": the grid's own vertices are up to ",
                                dvert, " off that shape and are LEFT THERE (moving them would ",
                                "break conformity with the neighbouring elements). The curved ",
                                "boundary interpolates them, so it is that far off the shape too.");
                         msg_rank = rank, suppress = mesh.msg_suppress)
        end
    end

    ncurved_tags == 0 && return nothing

    #-----------------------------------------------------------------------------
    # Linear blending weights. `lgl.ξ` runs from -1 to +1, so for a curved edge
    # sitting at index 1 of a direction the weight is 1 there and decays
    # linearly to 0 at index ngl, and vice versa. Written as a ratio of the end
    # values rather than (1∓ξ)/2 so it stays exact whatever the endpoint
    # convention of the interpolation nodes is.
    #-----------------------------------------------------------------------------
    ξ    = TFloat.(collect(lgl.ξ[1:ngl]))
    span = ξ[ngl] - ξ[1]
    wlo  = TFloat[(ξ[ngl] - ξ[k])/span for k = 1:ngl]      # 1 at index 1, 0 at ngl
    whi  = TFloat[(ξ[k] - ξ[1])/span   for k = 1:ngl]      # 0 at index 1, 1 at ngl

    #-----------------------------------------------------------------------------
    # Gordon-Hall blend into the interiors.
    #
    # A direction line of the element is on the curved boundary exactly when all
    # ngl of its nodes are. It cannot be a false positive: for ngl >= 3 the
    # interior nodes of an edge belong to that edge alone, so they are marked
    # only if that very edge was snapped.
    #-----------------------------------------------------------------------------
    nelem_curved = 0
    dmax_int     = 0.0
    for iel = 1:mesh.nelem

        jlo = true; jhi = true; ilo = true; ihi = true
        for k = 1:ngl
            jlo &= on_curved[mesh.connijk[iel, k, 1]]
            jhi &= on_curved[mesh.connijk[iel, k, ngl]]
            ilo &= on_curved[mesh.connijk[iel, 1, k]]
            ihi &= on_curved[mesh.connijk[iel, ngl, k]]
        end
        (jlo || jhi || ilo || ihi) || continue
        nelem_curved += 1

        for j = 2:ngl-1, i = 2:ngl-1
            ip = mesh.connijk[iel, i, j]
            dx = zero(TFloat); dy = zero(TFloat)
            if jlo
                p = mesh.connijk[iel, i, 1];   dx += wlo[j]*δx[p]; dy += wlo[j]*δy[p]
            end
            if jhi
                p = mesh.connijk[iel, i, ngl]; dx += whi[j]*δx[p]; dy += whi[j]*δy[p]
            end
            if ilo
                p = mesh.connijk[iel, 1, j];   dx += wlo[i]*δx[p]; dy += wlo[i]*δy[p]
            end
            if ihi
                p = mesh.connijk[iel, ngl, j]; dx += whi[i]*δx[p]; dy += whi[i]*δy[p]
            end
            δx[ip] = dx; δy[ip] = dy
            dmax_int = max(dmax_int, hypot(dx, dy))
        end
    end

    #-----------------------------------------------------------------------------
    # Commit. One pass over everything: nodes that were never touched carry
    # δ = 0 and so are bit-for-bit unchanged.
    #-----------------------------------------------------------------------------
    @inbounds for ip = 1:mesh.npoin
        mesh.x[ip] += δx[ip]
        mesh.y[ip] += δy[ip]
    end

    println_rank(string(" #   ", MPI.Allreduce(nelem_curved, MPI.SUM, comm),
                        " elements curved ; largest interior node correction = ",
                        MPI.Allreduce(dmax_int, MPI.MAX, comm));
                 msg_rank = rank, suppress = mesh.msg_suppress)

    #-----------------------------------------------------------------------------
    # Validity. Curving a boundary on a grid that is too coarse for it — the
    # element radially thinner than the sagitta of its own boundary arc — folds
    # the element and produces a negative Jacobian. Catch it here, where the
    # message can name the geometry, rather than letting it surface as a NaN
    # thousands of time steps later.
    #-----------------------------------------------------------------------------
    _check_curved_elements(mesh, ξ, on_curved, comm, rank)

    # The outer boundary may itself have been curved, so the bounding box that
    # mod_mesh_read_gmsh! computed from the straight-sided grid can be stale.
    mesh.xmax = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
    mesh.xmin = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
    mesh.ymax = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm)
    mesh.ymin = MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)

    return nothing
end


#---------------------------------------------------------------------------------
# Jacobian sign check on the elements that were touched, using the same
# cross-product form (Kopriva eq. 38) that build_metric_terms! will use:
#
#   J = X_ξ Y_η - Y_ξ X_η ,  differentiated at the nodes.
#
# Only the sign matters here, and it must agree with the sign the grid already
# had — an element that flipped is a folded element.
#---------------------------------------------------------------------------------
function _check_curved_elements(mesh::St_mesh, ξ::AbstractVector, on_curved, comm, rank)

    ngl = mesh.ngl
    D   = lagrange_nodal_derivative_matrix(ξ)

    nbad  = 0
    ratio = Inf                          # min |J| / max |J| over curved elements
    for iel = 1:mesh.nelem
        touched = false
        for k = 1:ngl
            touched |= on_curved[mesh.connijk[iel, k, 1]]   | on_curved[mesh.connijk[iel, k, ngl]] |
                       on_curved[mesh.connijk[iel, 1, k]]   | on_curved[mesh.connijk[iel, ngl, k]]
        end
        touched || continue

        Jmin = Inf; Jmax = 0.0; sgn = 0
        for j = 1:ngl, i = 1:ngl
            xξ = 0.0; yξ = 0.0; xη = 0.0; yη = 0.0
            for k = 1:ngl
                p = mesh.connijk[iel, k, j]
                xξ += D[i,k]*mesh.x[p]; yξ += D[i,k]*mesh.y[p]
                q = mesh.connijk[iel, i, k]
                xη += D[j,k]*mesh.x[q]; yη += D[j,k]*mesh.y[q]
            end
            J = xξ*yη - yξ*xη
            s = J > 0 ? 1 : (J < 0 ? -1 : 0)
            sgn == 0 && (sgn = s)
            if s != sgn                                      # J changed sign inside: folded
                nbad += 1
                break
            end
            Jmin = min(Jmin, abs(J)); Jmax = max(Jmax, abs(J))
        end
        Jmax > 0 && (ratio = min(ratio, Jmin/Jmax))
    end

    nbad = MPI.Allreduce(nbad, MPI.SUM, comm)
    if nbad > 0
        error(" # ERROR exact_geometry.jl: curving the boundary folded " * string(nbad) *
              " element(s): the Jacobian changes sign inside them.\n" *
              " #   The grid is too coarse normal to the curved boundary for its own\n" *
              " #   curvature. Refine the boundary layer, or drop the :exact_geometry entry.")
    end
    println_rank(string(" #   curved-element Jacobian check passed ; min|J|/max|J| within an element = ",
                        MPI.Allreduce(ratio, MPI.MIN, comm));
                 msg_rank = rank, suppress = mesh.msg_suppress)

    return nothing
end
