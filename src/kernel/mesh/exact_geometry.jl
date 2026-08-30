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
# Two things, and they live in two different places on purpose:
#
#   1. THE DECK names which boundaries are curved, and may carry a per-boundary
#      parameter for the case's own use:
#
#        :exact_geometry => ["circle_boundary"]                     # tag list
#        :exact_geometry => "circle_boundary"                       # one tag
#        :exact_geometry => Dict("circle_boundary" => (:circle, 1.0, 0.0, 0.2))
#
#      where <boundary tag> is the gmsh `Physical Curve` name. In the Dict form
#      the value is a SPEC: this file never looks inside it, it just hands it
#      back to the case. In the list form the spec is `nothing`.
#
#   2. THE CASE states the geometry, analytically, in
#      problems/<eqs>/<case>/user_exactGeo.jl:
#
#        # REQUIRED. Where does the point (x, y) belong on boundary `tag`?
#        # Return it unchanged for a tag this case does not handle.
#        user_exactGeo(tag, spec, x, y) -> (x, y)
#
#        # OPTIONAL. Called ONCE per tag, before any node is moved, with the
#        # LINEAR (gmsh vertex) point cloud of that boundary. Return the spec to
#        # hand to user_exactGeo — possibly a different one, fitted from the
#        # grid — or `nothing` to leave this boundary alone.
#        user_exactGeo_setup(tag, spec, xs, ys) -> spec or nothing
#
#      problems/CompEuler/shock_circle/user_exactGeo.jl is the worked example: a
#      circle, both stated outright and least-squares fitted from the grid (via
#      `je_fit_circle` below). Copy it and change the shape. Nothing about a
#      circle — or about any other particular geometry — is hardcoded here; this
#      file only knows how to ASK, blend and check.
#
# No :exact_geometry key, or an empty one, means "do nothing" — every existing
# case is unaffected. Naming a boundary without shipping a user_exactGeo.jl is
# an error, not a silent no-op.
#
# LIMITS
# ------
#   * 2D only. The 3D construction blends six FACE surfaces rather than four
#     edge curves, and Theorem 5 of the same paper shows the cross-product
#     metrics do NOT satisfy the metric identities there whatever the mapping
#     is — 3D needs the conservative curl form as well, so it is its own change.
#     NSD_1D and NSD_3D are explicit no-ops rather than silent ones. When 3D
#     lands, the hook grows a third coordinate; the 2D one stays as it is.
#
#   * ONE ELEMENT LAYER. Only the elements that own a curved edge are curved,
#     which is the standard construction (it is the mesh of Kopriva's Fig. 2)
#     and is what keeps the correction local and the grid conforming. The grid
#     must therefore already resolve the curvature: an element thinner, normal
#     to the wall, than the sagitta of its own boundary arc folds, and
#     _check_curved_elements stops the run rather than hand a negative Jacobian
#     to build_metric_terms!.
#
#   * REFINEMENT NEEDS THE GEOMETRY STATED, not fitted. Refinement happens on the
#     straight-sided Gridap model and this routine then runs on the result, which
#     is the right order — but it means each new boundary vertex arrives at the
#     MIDPOINT OF A CHORD, a sagitta inside the true curve (3.8e-3 on shock_circle
#     at :init_refine_lvl => 1). The snap handles that: vertices are moved too,
#     and the wall comes out at 3e-16 from the circle, refined or not. What CANNOT
#     survive it is a spec the case FITS from the grid, because a fit from
#     vertices that are half on the curve and half a sagitta inside it lands
#     exactly between the two clusters — shock_circle's user_exactGeo_setup
#     refuses it. On a refined grid state the geometry outright. Verified through
#     mod_mesh_mesh_driver at :init_refine_lvl => 1: 32 curved elements, wall at
#     3.3e-16, no fold. :ladapt (adapting DURING the run) is still untested.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# THE USER HOOKS.
#
# The geometry itself is the ONE thing this file does not know and must not
# guess, because it is different for every case: a circle here, an ellipse, an
# aerofoil, a bump, a bell, a spline through surveyed points. So it is asked
# for, from the case directory, exactly the way user_flux!/user_source!/
# user_analytic_solution already are:
#
#     problems/<eqs>/<case>/user_exactGeo.jl
#
# The definitions below are the DEFAULTS, in force for a case that ships no such
# file. `user_exactGeo_setup` passes the deck's spec straight through, and
# `user_exactGeo` refuses: naming a boundary in :exact_geometry and then not
# saying what shape it is is a mistake in the deck, and a silent no-op would
# leave the run quietly on the polygon it was trying to get off.
#
# A case's own user_exactGeo.jl defines these at the SAME signatures, so
# including it replaces them (see the note in src/run.jl on why the include is
# the whole mechanism). snap_nodes_to_exact_geometry! only calls them when
# inputs[:_has_exactgeo] says the case really shipped the file, so a definition
# left behind by a previously-run case in the same session can never be used.
#---------------------------------------------------------------------------------

#
# user_exactGeo(tag, spec, x, y) -> (x, y)
#
# THE analytical definition: given a point at or near boundary `tag`, return the
# corresponding point ON the exact geometry — for a closed shape, its nearest
# point. `spec` is whatever the deck wrote for this tag (`nothing` in the list
# form), after user_exactGeo_setup has had a look at it.
#
# Return (x, y) unchanged to leave a node where it is; that is also how a case
# declines a tag it does not handle.
#
# Theorem 3 does not care WHICH parametrization of the arc the N+1 nodes end up
# at — only that they lie on the curve, so that the interpolant through them is
# the P^N boundary. A nearest-point projection is the simple choice and the one
# the shipped example makes; anything onto-the-curve is admissible.
#
function user_exactGeo(tag, spec, x, y)
    error(" # ERROR exact_geometry.jl: :exact_geometry names the boundary \"" *
          String(tag) * "\", but this case does not say what shape it is.\n" *
          " #   Add problems/<eqs>/<case>/user_exactGeo.jl defining\n" *
          " #       user_exactGeo(tag, spec, x, y) -> (x, y)\n" *
          " #   (see problems/CompEuler/shock_circle/user_exactGeo.jl), or drop\n" *
          " #   the :exact_geometry entry.")
end

#
# user_exactGeo_setup(tag, spec, xs, ys) -> spec, or nothing to skip the tag
#
# Optional. Called once per boundary tag, before any node is moved, with the
# LINEAR (gmsh vertex) coordinates of that boundary gathered over all ranks —
# the geometry the mesh file itself asserts. Use it to fit the shape parameters
# from the grid (`je_fit_circle` below is there for exactly that), to validate
# what the deck asked for, or to say "not this one" by returning `nothing`.
#
# EVERY RANK CALLS IT, for every tag, in the same order — which is why a
# collective is allowed here and nowhere else, and why its verdict has to be the
# same on all of them: a `nothing` on some ranks only would skip the tag there
# and leave the others inside the next collective, hanging the run.
#
# The default hands the deck's spec back untouched.
#
user_exactGeo_setup(tag, spec, xs, ys) = spec


#---------------------------------------------------------------------------------
# je_fit_circle(xs, ys) -> (xc, yc, r, relative residual), or nothing
#
# Least-squares circle through a point cloud (Kåsa / algebraic fit). A NUMERICAL
# UTILITY, not a geometry: it is here rather than in a user file because it is
# an MPI collective, and a user_exactGeo_setup that wants the circle the grid
# itself states should call this instead of writing the reduction by hand. See
# problems/CompEuler/shock_circle/user_exactGeo.jl for the call.
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
# boundary. EVERY RANK MUST CALL IT, in the same order, for the same reason.
#
# The fourth return value is the largest radial residual relative to r: how far
# the fitted circle actually is from the vertices it was fitted to. A boundary
# that is not a circle shows up there, and the caller is expected to look — a
# silent bad fit would deform the wall onto a shape nobody asked for.
#
# Returns `nothing` when the cloud is degenerate (fewer than 3 points, or
# collinear — a straight boundary, for which the normal system is singular).
#
#---------------------------------------------------------------------------------
function je_fit_circle(xs::AbstractVector, ys::AbstractVector, comm = get_mpi_comm())

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

    return (TFloat(xc), TFloat(yc), TFloat(r), resid/r)
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
# Normalise whatever the deck wrote for :exact_geometry into a sorted
# [(tag::String, spec), ...].
#
# Accepted, all naming the same thing:
#
#   "circle_boundary"                                   one tag,  spec = nothing
#   :circle_boundary                                    ditto (Symbol keys are
#                                                       how a deck often spells
#                                                       a gmsh group)
#   ["wall", "hub"]  / ("wall", "hub")  / Set(...)      tags,     spec = nothing
#   Dict("wall" => (:circle, 1.0, 0.0, 0.2))            tag => spec
#
# SORTED BY TAG on purpose: user_exactGeo_setup may run an MPI collective, so
# every rank has to walk the boundaries in the same order, and a Dict's
# iteration order is not that.
#---------------------------------------------------------------------------------
_exact_geometry_entries(::Nothing) = Tuple{String,Any}[]

_exact_geometry_entries(spec::Union{AbstractString,Symbol}) =
    Tuple{String,Any}[(String(spec), nothing)]

_exact_geometry_entries(spec::AbstractDict) =
    sort!(Tuple{String,Any}[(String(k), v) for (k, v) in spec]; by = first)

_exact_geometry_entries(spec) =                       # any iterable of tags
    sort!(Tuple{String,Any}[(String(t), nothing) for t in spec]; by = first)


#---------------------------------------------------------------------------------
# Characteristic length of a boundary, from its own vertex cloud: the diagonal
# of its bounding box, reduced over the ranks that hold a piece of it. Used only
# to turn "this vertex has not moved" into a relative test.
#---------------------------------------------------------------------------------
function _boundary_scale(xs, ys, comm)
    lo = TFloat[isempty(xs) ?  Inf : minimum(xs), isempty(ys) ?  Inf : minimum(ys)]
    hi = TFloat[isempty(xs) ? -Inf : maximum(xs), isempty(ys) ? -Inf : maximum(ys)]
    lo = MPI.Allreduce(lo, MPI.MIN, comm)
    hi = MPI.Allreduce(hi, MPI.MAX, comm)
    d  = hypot(hi[1] - lo[1], hi[2] - lo[2])
    return isfinite(d) && d > 0 ? d : one(TFloat)
end


#---------------------------------------------------------------------------------
# Resolve one <boundary tag> => <shape spec> entry, by asking the case.
#
# All this does is call the case's user_exactGeo_setup with the boundary's own
# linear vertices and let it say what the spec really is — fitted, validated, or
# passed straight through — or decline the tag by returning `nothing`. It is a
# named function rather than an inline call so the caller reads the same as it
# did before, and so a case that declines a boundary is reported once, here.
#
# EVERY RANK CALLS THIS FOR EVERY TAG, in the same order: user_exactGeo_setup
# may (and shock_circle's does) run an MPI collective.
#---------------------------------------------------------------------------------
function _resolve_shape(tag::AbstractString, spec, xs, ys, comm, rank, suppress)

    shape = user_exactGeo_setup(tag, spec, xs, ys)

    if shape === nothing
        println_rank(string(" #   :exact_geometry \"", tag,
                            "\" skipped: user_exactGeo_setup declined it. Nothing was curved ",
                            "on that boundary.");
                     msg_rank = rank, suppress = suppress)
    end
    return shape
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

    entries = _exact_geometry_entries(get(inputs, :exact_geometry, nothing))
    isempty(entries) && return nothing

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    # The case must actually ship a user_exactGeo.jl. Checked through the deck
    # flag rather than a bare isdefined()/hasmethod(): user_exactGeo is a method
    # on the Jexpresso module, so once ANY case that defines it has run in this
    # session the definition outlives a switch to another case, and a case that
    # forgot the file would silently be curved onto its predecessor's geometry.
    # src/run.jl sets this from isfile(<case dir>/user_exactGeo.jl).
    if !get(inputs, :_has_exactgeo, false)
        error(" # ERROR exact_geometry.jl: :exact_geometry names " *
              string(length(entries)) * " boundary(ies) — " *
              join(("\"" * t * "\"" for (t, _) in entries), ", ") *
              " — but this case ships no user_exactGeo.jl,\n" *
              " #   so nothing says what shape they are. Add\n" *
              " #       problems/<eqs>/<case>/user_exactGeo.jl\n" *
              " #   defining  user_exactGeo(tag, spec, x, y) -> (x, y)  (copy\n" *
              " #   problems/CompEuler/shock_circle/user_exactGeo.jl), or drop the\n" *
              " #   :exact_geometry entry from user_inputs.jl.")
    end

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
    δx = zeros(TFloat, mesh.npoin)
    δy = zeros(TFloat, mesh.npoin)
    on_curved = falses(mesh.npoin)

    # `entries` is sorted by tag, so every rank walks them in the same order:
    # user_exactGeo_setup may run an MPI collective (a fit over the boundary
    # vertices) and must be entered in lockstep.
    ncurved_tags = 0
    for (tagstr, shape_spec) in entries

        edges = findall(t -> t !== nothing && String(t) == tagstr, mesh.bdy_edge_type)

        # Vertices of this boundary, for user_exactGeo_setup. Only the LINEAR
        # ones: they are the geometry the mesh file asserts; the high-order ones
        # are the chord interpolant we are about to correct.
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

        # "Already on the geometry" tolerance for a linear vertex, relative to
        # the size of the boundary itself. Measured from the vertex cloud rather
        # than asked of the shape, because the shape is the case's own object and
        # this file makes no assumption about what it is.
        vtol = 1.0e-10*_boundary_scale(xs, ys, comm)

        dmax = 0.0
        dvert = 0.0
        for iedge in edges
            for igl = 1:ngl
                ip = mesh.poin_in_bdy_edge[iedge, igl]
                on_curved[ip] = true
                xs_new, ys_new = user_exactGeo(tagstr, shape, mesh.x[ip], mesh.y[ip])
                dxi = xs_new - mesh.x[ip]
                dyi = ys_new - mesh.y[ip]
                d   = hypot(dxi, dyi)
                if igl == 1 || igl == ngl
                    d <= vtol && continue                  # already there
                    dvert = max(dvert, d)
                end
                δx[ip] = dxi; δy[ip] = dyi
                dmax = max(dmax, d)
            end
        end

        println_rank(string(" #   \"", tagstr, "\" -> user_exactGeo(", repr(shape), ")",
                            " ; ", nedges_glob, " edges ; largest node correction = ",
                            MPI.Allreduce(dmax, MPI.MAX, comm),
                            " (of which at a linear vertex: ",
                            MPI.Allreduce(dvert, MPI.MAX, comm), ")");
                     msg_rank = rank, suppress = mesh.msg_suppress)
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
    # FULL Gordon-Hall blend, corner terms included:
    #
    #   δ(ξ,η) = wlo(η)Γ₁(ξ) + whi(η)Γ₂(ξ) + wlo(ξ)Γ₃(η) + whi(ξ)Γ₄(η)
    #            - [ wlo(ξ)wlo(η)δ₁₁ + whi(ξ)wlo(η)δ₂₁
    #              + wlo(ξ)whi(η)δ₁₂ + whi(ξ)whi(η)δ₂₂ ]
    #
    # where Γ is the displacement along an element edge and δ.. the four corner
    # displacements. WHY THIS IS CONFORMING once the vertices move:
    #
    #   * a corner is one global node, so both elements sharing it use the same
    #     δ.. ;
    #   * on a CURVED boundary edge, Γ is the snap displacement — and a boundary
    #     edge belongs to exactly one element, so there is nobody to disagree
    #     with;
    #   * on ANY OTHER edge, Γ is defined below as the linear interpolation of
    #     that edge's own two corner displacements. It depends on nothing but
    #     those two global values, so the two elements sharing the edge compute
    #     it identically;
    #   * and the formula reproduces Γ exactly on each edge and δ.. exactly at
    #     each corner (the corner bracket cancels the edge terms there), so the
    #     three bullets above are the whole story: writing it at every one of
    #     the ngl² nodes is single-valued, from whichever element you look.
    #
    # A direction line is a curved boundary edge exactly when all ngl of its
    # nodes are flagged. That cannot be a false positive: for ngl >= 3 the
    # interior nodes of an edge belong to that edge alone, so they carry the
    # flag only if that very edge was snapped.
    #-----------------------------------------------------------------------------
    δx_gh = zeros(TFloat, mesh.npoin)          # the blended field, built fresh:
    δy_gh = zeros(TFloat, mesh.npoin)          # δx/δy above stay the BOUNDARY data
    written = falses(mesh.npoin)               # conformity guard, see below
    Γx    = Array{TFloat}(undef, 4, ngl)       # per-element edge displacements,
    Γy    = Array{TFloat}(undef, 4, ngl)       # order: j=1, j=ngl, i=1, i=ngl
    mismatch = 0.0                             # worst disagreement between two
                                               # elements over a shared node

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

        # The four corners, and their displacements.
        c11 = mesh.connijk[iel, 1, 1];     c21 = mesh.connijk[iel, ngl, 1]
        c12 = mesh.connijk[iel, 1, ngl];   c22 = mesh.connijk[iel, ngl, ngl]
        d11x = δx[c11]; d11y = δy[c11];    d21x = δx[c21]; d21y = δy[c21]
        d12x = δx[c12]; d12y = δy[c12];    d22x = δx[c22]; d22y = δy[c22]

        # Nothing to do unless this element owns a curved edge or has a moved
        # corner. (A moved corner reaches elements that touch the boundary only
        # at a point — they have to be blended too, or their edge into that
        # corner would tear away from it.)
        anycorner = (d11x != 0) | (d11y != 0) | (d21x != 0) | (d21y != 0) |
                    (d12x != 0) | (d12y != 0) | (d22x != 0) | (d22y != 0)
        (jlo || jhi || ilo || ihi || anycorner) || continue
        nelem_curved += 1

        # Edge displacements Γ. Curved boundary edge -> the snapped δ; any other
        # edge -> linear in its own two corner δ, which is what both of its
        # elements compute.
        for k = 1:ngl
            if jlo
                p = mesh.connijk[iel, k, 1];   Γx[1,k] = δx[p];  Γy[1,k] = δy[p]
            else
                Γx[1,k] = wlo[k]*d11x + whi[k]*d21x
                Γy[1,k] = wlo[k]*d11y + whi[k]*d21y
            end
            if jhi
                p = mesh.connijk[iel, k, ngl]; Γx[2,k] = δx[p];  Γy[2,k] = δy[p]
            else
                Γx[2,k] = wlo[k]*d12x + whi[k]*d22x
                Γy[2,k] = wlo[k]*d12y + whi[k]*d22y
            end
            if ilo
                p = mesh.connijk[iel, 1, k];   Γx[3,k] = δx[p];  Γy[3,k] = δy[p]
            else
                Γx[3,k] = wlo[k]*d11x + whi[k]*d12x
                Γy[3,k] = wlo[k]*d11y + whi[k]*d12y
            end
            if ihi
                p = mesh.connijk[iel, ngl, k]; Γx[4,k] = δx[p];  Γy[4,k] = δy[p]
            else
                Γx[4,k] = wlo[k]*d21x + whi[k]*d22x
                Γy[4,k] = wlo[k]*d21y + whi[k]*d22y
            end
        end

        for j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j]
            dx = wlo[j]*Γx[1,i] + whi[j]*Γx[2,i] + wlo[i]*Γx[3,j] + whi[i]*Γx[4,j] -
                 (wlo[i]*wlo[j]*d11x + whi[i]*wlo[j]*d21x +
                  wlo[i]*whi[j]*d12x + whi[i]*whi[j]*d22x)
            dy = wlo[j]*Γy[1,i] + whi[j]*Γy[2,i] + wlo[i]*Γy[3,j] + whi[i]*Γy[4,j] -
                 (wlo[i]*wlo[j]*d11y + whi[i]*wlo[j]*d21y +
                  wlo[i]*whi[j]*d12y + whi[i]*whi[j]*d22y)
            # CONFORMITY GUARD. A node on a shared edge is blended by both of
            # its elements, and the argument above says they must agree. If they
            # ever do not, the grid tears at that edge and every symptom of it
            # appears far downstream (a DSS summing values at two different
            # points), so measure the disagreement here rather than trust the
            # argument. Costs one comparison per element-node at grid-build time.
            if written[ip]
                mismatch = max(mismatch, hypot(dx - δx_gh[ip], dy - δy_gh[ip]))
            end
            written[ip] = true
            δx_gh[ip] = dx; δy_gh[ip] = dy
            (1 < i < ngl && 1 < j < ngl) && (dmax_int = max(dmax_int, hypot(dx, dy)))
        end
    end

    # Round-off is the only disagreement allowed; anything larger is a torn grid.
    mismatch_g = MPI.Allreduce(mismatch, MPI.MAX, comm)
    scale_g    = MPI.Allreduce(maximum(hypot.(δx_gh, δy_gh); init = 0.0), MPI.MAX, comm)
    if mismatch_g > max(1.0e-10*scale_g, 1.0e-13)
        error(" # ERROR exact_geometry.jl: the curved grid is NOT conforming — two elements\n" *
              " #   disagree by " * string(mismatch_g) * " about where a node they share must go\n" *
              " #   (largest displacement in the blend: " * string(scale_g) * ").\n" *
              " #   This is a bug in the Gordon-Hall blend, not in the case: please report it.")
    end

    #-----------------------------------------------------------------------------
    # Commit. One pass over everything: nodes no element blended carry δ = 0 and
    # so are bit-for-bit unchanged.
    #-----------------------------------------------------------------------------
    @inbounds for ip = 1:mesh.npoin
        mesh.x[ip] += δx_gh[ip]
        mesh.y[ip] += δy_gh[ip]
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
