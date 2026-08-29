#---------------------------------------------------------------------------------
# test/test_exact_geometry.jl — 2D curved-boundary snapping
# (src/kernel/mesh/exact_geometry.jl).
#
#   julia test/test_exact_geometry.jl
#
# NO `using Jexpresso`. Like test_cubed_sphere_maps.jl, this includes the source
# directly; the handful of names exact_geometry.jl borrows from the package
# (TFloat, MPI, println_rank, St_mesh, NSD_*) are stubbed below, so the numerics
# stay checkable without building a 130-package environment. The stub St_mesh
# carries exactly the fields the routine touches — if it grows a dependency on
# another field, this file stops compiling, which is the intended alarm.
#
# The grid under test is an annulus of nr x nt quads around a circle of radius
# R: the shock_circle topology in miniature, built STRAIGHT-SIDED the way
# mod_mesh_read_gmsh! builds one (linear vertices on the circle, LGL nodes
# interpolated on the straight quad), so the "before" state is exactly what the
# solver would otherwise see, and the exact answer is known in closed form.
#
# What is being pinned down, in the order it matters:
#
#   1. THE GRID STAYS CONFORMING. This is the property the whole construction is
#      arranged around and the one whose failure would be silent and fatal: a
#      node is shared between elements, so if the curving moved a node that a
#      NON-curved neighbour also owns, the two elements would disagree about
#      where their shared edge is and the DSS would be summing values at
#      different points. It cannot happen here because the transfinite blend of
#      an edge displacement vanishes identically on the two edges adjacent to
#      that edge — which is only true because the element corners (the linear
#      gmsh vertices) are left exactly where they were. Asserted directly:
#      every node of every element that does not touch the circle is bit-for-bit
#      unmoved.
#
#   2. THE ELEMENTS DO NOT FOLD. Curving a boundary on a grid too coarse normal
#      to it turns the element inside out. The Jacobian is checked over the
#      whole grid with the same cross-product form build_metric_terms! uses.
#
#   3. FREE-STREAM PRESERVATION SURVIVES. Kopriva (2006) Theorem 3: in 2D the
#      cross-product metrics satisfy the discrete metric identities iff the
#      mapping is in P^N. Curving moves nodes, and the mapping is the degree-N
#      interpolant through them, so it stays in P^N — but that is an argument,
#      and this is the measurement. Form II of the DMIs, eq. (30), is evaluated
#      before and after and must sit at round-off both times.
#
#   4. THE BOUNDARY IS ACTUALLY THE CIRCLE NOW. Not just at the nodes (that is
#      trivially true after a projection) but BETWEEN them: the degree-N
#      interpolant through the snapped edge nodes is sampled densely and its
#      distance from the true circle measured. Straight-sided converges at
#      O(h²) and stops there however large :nop is — that is the defect being
#      fixed. Isoparametric degree-N edges converge at O(h^{N+1}) or better.
#
#   5. THE CURVED-EDGE SURFACE METRICS. build_metric_terms!(…, NSD_2D) gets Jef
#      and the wall normal from the edge tangent of the same interpolant. On a
#      straight edge that must reproduce, to round-off, the chord/2 and the
#      two-point normal it replaces (this is what makes the change a no-op for
#      every existing case); on the curved edge it must give the exact circle
#      normal, which the two-point difference does not.
#
#   6. The circle fit accepts a circle and REFUSES a straight boundary. The fit
#      is what runs when a deck writes `:circle` without coordinates, and a
#      silent bad fit would deform the grid onto a shape nobody asked for.
#---------------------------------------------------------------------------------

using Test
using LinearAlgebra
using Printf

#---------------------------------------------------------------------------------
# Package stubs — the whole of what exact_geometry.jl borrows from Jexpresso.
#---------------------------------------------------------------------------------
const TFloat = Float64
const TInt   = Int64

module MPI                                  # serial stand-in: Allreduce is identity
    struct Op; f::Function; end
    const SUM = Op(+); const MIN = Op(min); const MAX = Op(max)
    Comm_rank(c) = 0
    Allreduce(x, ::Op, c) = x
end
get_mpi_comm() = nothing
println_rank(args...; msg_rank = 0, suppress = false) = nothing

struct NSD_1D end
struct NSD_2D end
struct NSD_3D end

mutable struct St_mesh
    nsd::Int; nop::Int; ngl::Int
    npoin::Int; nelem::Int; nedges_bdy::Int
    x::Vector{Float64}; y::Vector{Float64}
    connijk::Array{Int,3}
    poin_in_bdy_edge::Array{Int,2}
    bdy_edge_type::Vector{Union{Nothing,String}}
    bdy_edge_in_elem::Vector{Int}
    SD::Any
    msg_suppress::Bool
    xmin::Float64; xmax::Float64; ymin::Float64; ymax::Float64
end

include(joinpath(@__DIR__, "..", "src", "kernel", "mesh", "exact_geometry.jl"))


#---------------------------------------------------------------------------------
# Legendre-Gauss-Lobatto nodes on [-1,1]: Newton on (1-x²)P'_N, started from the
# Chebyshev-Lobatto points. Same nodes basis_structs_ξ_ω! produces for "lgl".
#---------------------------------------------------------------------------------
function lgl_nodes(N::Int)
    N == 1 && return [-1.0, 1.0]
    x    = [-cos(pi*i/N) for i = 0:N]
    P    = zeros(N+1, N+1)
    xold = 2 .* x
    while maximum(abs.(x .- xold)) > 1.0e-15
        xold = copy(x)
        P[:,1] .= 1.0; P[:,2] .= x
        for k = 2:N
            @. P[:,k+1] = ((2k-1)*x*P[:,k] - (k-1)*P[:,k-1])/k
        end
        @. x = xold - (x*P[:,N+1] - P[:,N])/((N+1)*P[:,N+1])
    end
    x[1] = -1.0; x[end] = 1.0
    return x
end

# LGL quadrature weights, w_k = 2/(N(N+1) P_N(ξ_k)²).
function lgl_weights(N::Int)
    ξ = lgl_nodes(N)
    P = ones(N+1); Pm = zeros(N+1)
    for k = 1:N
        P, Pm = ((2k-1) .* ξ .* P .- (k-1) .* Pm) ./ k, P
    end
    return 2.0 ./ (N*(N+1) .* P.^2)
end

struct LGL; ξ::Vector{Float64}; end

const R, XC, YC = 0.2, 1.0, 0.0             # the circle every test below uses

#---------------------------------------------------------------------------------
# Straight-sided annulus of nr x nt quads, R <= r <= Rout, inner ring of edges
# tagged "circle". Index i runs radially (i == 1 on the circle), j azimuthally.
#---------------------------------------------------------------------------------
function build_annulus(N, nr, nt; Rout = 0.6)
    ngl = N + 1
    ξ   = lgl_nodes(N)
    id  = Dict{NTuple{2,Int},Int}()
    xs  = Float64[]; ys = Float64[]
    function push_node!(x, y)                    # de-duplicate shared nodes
        k = (round(Int, x/1.0e-12), round(Int, y/1.0e-12))
        haskey(id, k) && return id[k]
        push!(xs, x); push!(ys, y); id[k] = length(xs); return length(xs)
    end

    connijk = zeros(Int, nr*nt, ngl, ngl)
    bdy_p   = Int[]; bdy_e = Int[]
    rv(a) = R + (Rout - R)*a/nr
    θv(b) = 2pi*b/nt

    iel = 0
    for a = 0:nr-1, b = 0:nt-1
        iel += 1
        V = ntuple(4) do c
            (ra, tb) = ((a, b), (a+1, b), (a+1, b+1), (a, b+1))[c]
            (XC + rv(ra)*cos(θv(tb)), YC + rv(ra)*sin(θv(tb)))
        end
        for j = 1:ngl, i = 1:ngl                 # bilinear on the straight quad
            u = (ξ[i] + 1)/2; v = (ξ[j] + 1)/2
            x = (1-u)*(1-v)*V[1][1] + u*(1-v)*V[2][1] + u*v*V[3][1] + (1-u)*v*V[4][1]
            y = (1-u)*(1-v)*V[1][2] + u*(1-v)*V[2][2] + u*v*V[3][2] + (1-u)*v*V[4][2]
            connijk[iel, i, j] = push_node!(x, y)
        end
        if a == 0
            push!(bdy_e, iel)
            append!(bdy_p, [connijk[iel, 1, j] for j = 1:ngl])
        end
    end

    nbdy = length(bdy_e)
    mesh = St_mesh(2, N, ngl, length(xs), nr*nt, nbdy, xs, ys, connijk,
                   collect(reshape(bdy_p, ngl, nbdy)'),
                   Union{Nothing,String}["circle" for _ = 1:nbdy], bdy_e,
                   NSD_2D(), true,
                   minimum(xs), maximum(xs), minimum(ys), maximum(ys))
    return mesh, LGL(ξ)
end

curve!(mesh, lgl; shape = :circle) =
    snap_nodes_to_exact_geometry!(mesh, lgl,
        Dict{Symbol,Any}(:exact_geometry => Dict("circle" => shape)), mesh.SD)

# Metrics of one element exactly as build_metric_terms!(…, NSD_2D) forms them:
# by differentiating the degree-N nodal interpolant. Kopriva eq. (38).
function element_metrics(mesh, iel, D)
    ngl = mesh.ngl
    Ja1 = zeros(2, ngl, ngl); Ja2 = zeros(2, ngl, ngl); J = zeros(ngl, ngl)
    for j = 1:ngl, i = 1:ngl
        xξ = yξ = xη = yη = 0.0
        for k = 1:ngl
            p = mesh.connijk[iel, k, j]; xξ += D[i,k]*mesh.x[p]; yξ += D[i,k]*mesh.y[p]
            q = mesh.connijk[iel, i, k]; xη += D[j,k]*mesh.x[q]; yη += D[j,k]*mesh.y[q]
        end
        Ja1[1,i,j] =  yη; Ja1[2,i,j] = -xη
        Ja2[1,i,j] = -yξ; Ja2[2,i,j] =  xξ
        J[i,j] = xξ*yη - yξ*xη
    end
    return J, Ja1, Ja2
end

# Form II of the discrete metric identities, Kopriva eq. (30), normalised by the
# size of the metric terms so the number is dimensionless.
function dmi_residual(mesh, D)
    ngl = mesh.ngl; worst = 0.0; scale = 0.0
    for iel = 1:mesh.nelem
        _, Ja1, Ja2 = element_metrics(mesh, iel, D)
        scale = max(scale, maximum(abs, Ja1) + maximum(abs, Ja2))
        for n = 1:2, j = 1:ngl, i = 1:ngl
            s = 0.0
            for k = 1:ngl
                s += D[i,k]*Ja1[n,k,j] + D[j,k]*Ja2[n,i,k]
            end
            worst = max(worst, abs(s))
        end
    end
    return worst/scale
end

# Barycentric evaluation of the degree-N interpolant through (px, py) at ξ = s.
function eval_edge(ξ, px, py, s)
    n = length(ξ)
    for k = 1:n
        abs(s - ξ[k]) < 1.0e-14 && return (px[k], py[k])
    end
    X = Y = W = 0.0
    for k = 1:n
        w = 1.0
        for j = 1:n; j == k || (w /= (ξ[k] - ξ[j])); end
        w /= (s - ξ[k])
        X += w*px[k]; Y += w*py[k]; W += w
    end
    return (X/W, Y/W)
end

# Sup distance between the discrete boundary and the true circle.
function boundary_error(N, nt; curved, nr = 2)
    m, l = build_annulus(N, nr, nt)
    curved && curve!(m, l)
    worst = 0.0
    for ie = 1:m.nedges_bdy
        px = [m.x[m.poin_in_bdy_edge[ie,k]] for k = 1:m.ngl]
        py = [m.y[m.poin_in_bdy_edge[ie,k]] for k = 1:m.ngl]
        for s in range(-1, 1; length = 201)
            X, Y = eval_edge(l.ξ, px, py, s)
            worst = max(worst, abs(hypot(X - XC, Y - YC) - R))
        end
    end
    return worst
end


@testset verbose = true "exact geometry (2D curved boundaries)" begin

    @testset "nodal derivative matrix" begin
        for N in (2, 4, 7)
            ξ = lgl_nodes(N); D = lagrange_nodal_derivative_matrix(ξ)
            @test maximum(abs, D*ones(N+1))           < 1.0e-12         # d/dξ const
            @test maximum(abs, D*ξ .- 1.0)            < 1.0e-12         # d/dξ ξ
            # exact for every polynomial of degree <= N, and only those
            N >= 3 && @test maximum(abs, D*(ξ.^3) .- 3 .* ξ.^2) < 1.0e-11
        end
    end

    @testset "circle fit" begin
        m, _ = build_annulus(4, 2, 16)
        xs = Float64[]; ys = Float64[]
        for ie = 1:m.nedges_bdy, ig in (1, m.ngl)
            p = m.poin_in_bdy_edge[ie, ig]; push!(xs, m.x[p]); push!(ys, m.y[p])
        end
        fit = _fit_circle_mpi(xs, ys, nothing)
        @test fit !== nothing
        c, resid = fit
        @test c.xc ≈ XC atol = 1.0e-12
        @test c.yc ≈ YC atol = 1.0e-12
        @test c.r  ≈ R  atol = 1.0e-12
        @test resid < 1.0e-12

        # A straight boundary is not a circle: the fit must say so rather than
        # return a circle of enormous radius and deform the wall onto it.
        @test _fit_circle_mpi([0.0, 1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 0.0], nothing) === nothing
        @test _fit_circle_mpi([0.0, 1.0], [0.0, 1.0], nothing) === nothing      # too few

        # A square is not a circle either: the fit succeeds but the residual is
        # O(1) and _resolve_shape rejects it.
        sq = _fit_circle_mpi([0.0,1.0,1.0,0.0,0.5], [0.0,0.0,1.0,1.0,0.0], nothing)
        @test sq === nothing || sq[2] > 1.0e-6
        @test _resolve_shape("sq", :circle, [0.0,1.0,1.0,0.0,0.5], [0.0,0.0,1.0,1.0,0.0],
                             nothing, 0, true) === nothing
    end

    @testset "explicit shape spec" begin
        m, l = build_annulus(4, 2, 16)
        curve!(m, l; shape = (:circle, XC, YC, R))
        for ie = 1:m.nedges_bdy, ig = 1:m.ngl
            p = m.poin_in_bdy_edge[ie, ig]
            @test hypot(m.x[p] - XC, m.y[p] - YC) ≈ R atol = 1.0e-14
        end
        @test_throws ErrorException curve!(build_annulus(4, 2, 8)...; shape = :ellipse)

        # A deck may name the gmsh group with a Symbol just as well as a String.
        m2, l2 = build_annulus(4, 2, 16)
        snap_nodes_to_exact_geometry!(m2, l2,
            Dict{Symbol,Any}(:exact_geometry => Dict(:circle => (:circle, XC, YC, R))), m2.SD)
        @test m2.x == m.x && m2.y == m.y
    end

    @testset "no-op when nothing is asked for" begin
        m, l = build_annulus(4, 2, 16)
        x0, y0 = copy(m.x), copy(m.y)
        snap_nodes_to_exact_geometry!(m, l, Dict{Symbol,Any}(), m.SD)
        @test m.x == x0 && m.y == y0
        # a tag no boundary carries is reported, not applied
        snap_nodes_to_exact_geometry!(m, l,
            Dict{Symbol,Any}(:exact_geometry => Dict("typo" => :circle)), m.SD)
        @test m.x == x0 && m.y == y0
        # nop = 1: no high-order nodes exist, so there is nothing to move
        m1, l1 = build_annulus(1, 2, 16)
        z0, w0 = copy(m1.x), copy(m1.y)
        curve!(m1, l1)
        @test m1.x == z0 && m1.y == w0
    end

    @testset "snapped grid" begin
        N, nr, nt = 4, 4, 16
        mesh, lgl = build_annulus(N, nr, nt)
        D  = lagrange_nodal_derivative_matrix(lgl.ξ)
        x0, y0 = copy(mesh.x), copy(mesh.y)
        dmi0 = dmi_residual(mesh, D)
        curve!(mesh, lgl)

        @testset "boundary nodes land on the circle" begin
            off_before = maximum(abs(hypot(x0[mesh.poin_in_bdy_edge[ie,k]] - XC,
                                           y0[mesh.poin_in_bdy_edge[ie,k]] - YC) - R)
                                 for ie = 1:mesh.nedges_bdy, k = 1:mesh.ngl)
            off_after  = maximum(abs(hypot(mesh.x[mesh.poin_in_bdy_edge[ie,k]] - XC,
                                           mesh.y[mesh.poin_in_bdy_edge[ie,k]] - YC) - R)
                                 for ie = 1:mesh.nedges_bdy, k = 1:mesh.ngl)
            @test off_before > 1.0e-3            # the straight-sided defect is real
            @test off_after  < 1.0e-12
        end

        @testset "grid stays conforming" begin
            # (a) an element that does not touch the circle keeps every node.
            for iel = 1:mesh.nelem
                iel in mesh.bdy_edge_in_elem && continue
                for j = 1:mesh.ngl, i = 1:mesh.ngl
                    p = mesh.connijk[iel, i, j]
                    @test mesh.x[p] == x0[p]
                    @test mesh.y[p] == y0[p]
                end
            end
            # (b) the linear vertices — the element corners — never move, which
            #     is what kills the Gordon-Hall corner terms and makes (a) hold.
            for iel = 1:mesh.nelem, i in (1, mesh.ngl), j in (1, mesh.ngl)
                p = mesh.connijk[iel, i, j]
                @test mesh.x[p] == x0[p]
                @test mesh.y[p] == y0[p]
            end
            # (c) nothing moved outside the one element layer on the boundary.
            for p = 1:mesh.npoin
                (mesh.x[p] == x0[p] && mesh.y[p] == y0[p]) && continue
                @test hypot(x0[p] - XC, y0[p] - YC) <= R + (0.6 - R)/nr + 1.0e-12
            end
        end

        @testset "no element folds" begin
            for iel = 1:mesh.nelem
                J, _, _ = element_metrics(mesh, iel, D)
                @test minimum(J) > 0.0
            end
        end

        @testset "discrete metric identities (Kopriva Thm 3)" begin
            @test dmi0 < 1.0e-12
            @test dmi_residual(mesh, D) < 1.0e-12
        end
    end

    @testset "a fold is caught, not shipped" begin
        # One azimuthal element spanning 120° against a radial thickness far
        # smaller than the arc's sagitta: the blend has to turn the element
        # inside out, and the validity check must say so.
        @test_throws ErrorException begin
            m, l = build_annulus(4, 1, 3; Rout = R + 0.002)
            curve!(m, l)
        end
    end

    @testset "boundary geometry converges at the isoparametric rate" begin
        # Straight-sided is O(h²) whatever N is — the defect. Curved is at least
        # O(h^{N+1}); for a circle sampled radially it comes out one order
        # better still, so the assertion is >= N+1 with a little slack.
        for N in (3, 4)
            es = [boundary_error(N, nt; curved = false) for nt in (16, 32)]
            ec = [boundary_error(N, nt; curved = true)  for nt in (16, 32)]
            @test log2(es[1]/es[2]) ≈ 2.0 atol = 0.1          # stuck at 2nd order
            @test log2(ec[1]/ec[2]) >= N + 0.8                # isoparametric
            @test ec[2] < es[2]/1.0e2                         # and far smaller
        end
    end

    @testset "curved-edge surface metrics" begin
        # The Jef / normal formula that build_metric_terms!(…, NSD_2D) now uses.
        tangent_metrics(m, l, ie) = begin
            D = lagrange_nodal_derivative_matrix(l.ξ)
            ngl = m.ngl
            Jef = zeros(ngl); nx = zeros(ngl); ny = zeros(ngl)
            for k = 1:ngl
                tx = ty = 0.0
                for i = 1:ngl
                    p = m.poin_in_bdy_edge[ie, i]
                    tx += D[k,i]*m.x[p]; ty += D[k,i]*m.y[p]
                end
                Jef[k] = hypot(tx, ty); nx[k] = ty/Jef[k]; ny[k] = -tx/Jef[k]
            end
            return Jef, nx, ny
        end

        # (a) On a STRAIGHT edge it must reproduce the chord/2 and the two-point
        #     normal it replaces — this is what makes the change inert for every
        #     existing case.
        m, l = build_annulus(4, 2, 16)
        for ie = 1:m.nedges_bdy
            Jef, nx, ny = tangent_metrics(m, l, ie)
            p1 = m.poin_in_bdy_edge[ie, 1]; p2 = m.poin_in_bdy_edge[ie, m.ngl]
            chord = hypot(m.x[p1] - m.x[p2], m.y[p1] - m.y[p2])
            @test all(j -> isapprox(Jef[j], chord/2; atol = 1.0e-12), 1:m.ngl)
            dx = m.x[p2] - m.x[p1]; dy = m.y[p2] - m.y[p1]
            for k = 1:m.ngl
                @test nx[k] ≈  dy/chord atol = 1.0e-12
                @test ny[k] ≈ -dx/chord atol = 1.0e-12
            end
        end

        # (b) On the CURVED edge the normal must be the exact circle normal, and
        #     Jef the exact arc |dX/dξ| = R·dθ/dξ. The two-point difference that
        #     this replaces is only O(h) accurate, which is asserted too — it is
        #     the reason the surface metrics had to change with the grid.
        mc, lc = build_annulus(4, 2, 16)
        curve!(mc, lc)
        worst_exact = 0.0; worst_2pt = 0.0
        for ie = 1:mc.nedges_bdy
            Jef, nx, ny = tangent_metrics(mc, lc, ie)
            for k = 1:mc.ngl
                p  = mc.poin_in_bdy_edge[ie, k]
                ex = (mc.x[p] - XC)/R; ey = (mc.y[p] - YC)/R      # exact unit normal
                worst_exact = max(worst_exact, min(hypot(nx[k]-ex, ny[k]-ey),
                                                   hypot(nx[k]+ex, ny[k]+ey)))
                q  = mc.poin_in_bdy_edge[ie, k < mc.ngl ? k+1 : k-1]
                dx = mc.x[p] - mc.x[q]; dy = mc.y[p] - mc.y[q]; d = hypot(dx, dy)
                worst_2pt = max(worst_2pt, min(hypot(dy/d - ex, -dx/d - ey),
                                               hypot(dy/d + ex, -dx/d + ey)))
            end
            # Jef is |dX/dξ|, and radial projection does NOT make θ linear in
            # ξ, so it varies along the arc by a few percent — which is the
            # whole point: the chord/2 it replaces was constant. What must be
            # exact is the quantity Jef exists to produce, the arc length
            # ∫|dX/dξ|dξ collocated on the LGL nodes, against R·Δθ.
            p1 = mc.poin_in_bdy_edge[ie, 1]; p2 = mc.poin_in_bdy_edge[ie, mc.ngl]
            θ1 = atan(mc.y[p1] - YC, mc.x[p1] - XC)
            θ2 = atan(mc.y[p2] - YC, mc.x[p2] - XC)
            Δθ = abs(rem(θ2 - θ1, 2pi, RoundNearest))
            @test sum(lgl_weights(mc.nop) .* Jef) ≈ R*Δθ rtol = 1.0e-6
            @test !all(j -> isapprox(Jef[j], Jef[1]; rtol = 1.0e-3), 1:mc.ngl)
            @test all(j -> isapprox(Jef[j], R*Δθ/2; rtol = 0.05), 1:mc.ngl)
        end
        @test worst_exact < 1.0e-4          # tangent of the isoparametric edge
        @test worst_2pt   > 1.0e-2          # the two-point normal it replaces
        @test worst_2pt   > 50*worst_exact
    end
end
