#---------------------------------------------------------------------------------
# sphere_metrics.jl
#
# METRIC TERMS OF THE 2-D MANIFOLD EMBEDDED IN 3-D, and the (diagonal) mass
# matrix that goes with them. This is the piece that turns the F, G, H of
# user_flux.jl into a SURFACE DIVERGENCE, i.e. what was missing before the
# equations could be advanced in time.
#
# WHY THE FLAT 2-D METRICS DO NOT DO
# ----------------------------------
# On a flat 2-D grid the map (ξ,η) → (x,y) is square and invertible, and
# metric_terms.jl stores its inverse as dξdx, dξdy, dηdx, dηdy. On a shell the
# map is (ξ,η) → (x,y,z): a 3×2 Jacobian, which has no inverse. What replaces
# it is the CONTRAVARIANT BASIS of the surface.
#
# THE CONSTRUCTION
# ----------------
# At every LGL node of every element, differentiate the nodal coordinates with
# the LGL differentiation matrix to get the covariant (tangent) vectors
#
#   a_ξ = ∂x/∂ξ ,   a_η = ∂x/∂η        (3-vectors, tangent to the shell)
#
# From those, two DIFFERENT forms of the contravariant basis are available, and
# both are implemented here — pick one with
#
#   :sphere_metrics => :curl_invariant   (default)
#   :sphere_metrics => :cross_product
#
# in the case's user_inputs.jl. They agree in the continuum and differ, at any
# finite polynomial order, by the interpolation error of the sphere. Which one
# is used matters, and Section 3.2 of Kelly, Alves, Eckermann et al. (JCP 552,
# 2026, 114683) compares exactly these two (their CP and CI) inside NUMA and
# NEPTUNE, concluding that the curl-invariant form is the one to prefer for
# global-scale spectral-element simulation.
#
# (1) CROSS-PRODUCT FORM  (:cross_product) — the textbook one, Eq. (21) of
#     Kelly et al. The cross product of the tangents is normal to the surface
#     and its length is the surface Jacobian, i.e. the area element:
#
#       n = a_ξ × a_η ,   J = |n| ,   n̂ = n/J
#
#     and the contravariant basis is
#
#       a¹ = (a_η × n)/J² ,   a² = (n × a_ξ)/J² .
#
# (2) CURL-INVARIANT FORM  (:curl_invariant) — Kopriva, "Metric identities and
#     the discontinuous spectral element method on curvilinear meshes",
#     J. Sci. Comput. 26(3), 301-327 (2006), Eq. (15); Eq. (23) of Kelly et al.
#     In three dimensions the metric terms are written as a CURL,
#
#       J aⁱ = ½ [ ∂_ξk (∂x/∂ξʲ × x) − ∂_ξj (∂x/∂ξᵏ × x) ] ,  (i,j,k) cyclic,
#
#     which is divergence-free by inspection, so the metric identities
#     Σᵢ ∂_ξi (J aⁱ) = 0 survive DISCRETELY: the discrete divergence of a
#     discrete curl vanishes identically whenever the same differentiation
#     matrix is used for both, whereas the cross-product form only satisfies
#     them to the order of the approximation. That is Kopriva's whole point and
#     it is what "constant-state preservation" rests on.
#
#     On a 2-D manifold there is no third reference coordinate to take that
#     curl in, so one has to be supplied. The shell has an obvious one: extend
#     the surface radially, x(ξ,η,ζ) = (ζ/R) x(ξ,η), so that ζ = r is the
#     radius and the shell is the level set ζ = R. Feeding that into the
#     formula above and evaluating on the shell collapses the ζ-derivatives
#     analytically and leaves, with x̂ = x/|x| the unit radial,
#
#       J a¹ = a_η × x̂ ,   J a² = x̂ × a_ξ ,   J = x̂·(a_ξ × a_η) ,
#
#       J a³ = J n̂ = (R/2) [ ∂_η (a_ξ × x̂) − ∂_ξ (a_η × x̂) ]      (*)
#
#     — the first two are pointwise products of nodal values, the third is the
#     surviving discrete curl. (*) is NOT a_ξ × a_η discretely, and that is the
#     entire content of the method: the three of them together satisfy the
#     metric identity of the extended map to machine precision,
#
#       ∂_ξ(J a¹) + ∂_η(J a²) + (2/R) (J n̂) = 0 ,                  (**)
#
#     which is the curved-surface statement — the 2/R is the mean curvature of
#     the sphere; see M5 of check_sphere_metrics for why a flat ∇·(J aⁱ) = 0
#     would be the wrong identity to assert here. J n̂ is stored alongside the
#     rest (Jnx, Jny, Jnz) precisely so that (**) can be checked, and so that
#     the RHS has the curl-consistent normal available rather than a
#     recomputed cross product.
#
#     Two by-products come for free and are worth naming, because the
#     cross-product form has neither: J a¹ and J a² are cross products with the
#     EXACT radial direction, so aⁱ·x̂ = 0 to round-off — the contravariant
#     basis is tangent to the true sphere, not merely to the polynomial
#     interpolant of it — and aⁱ·a_j = δⁱⱼ holds exactly by the choice
#     J = x̂·(a_ξ × a_η) that goes with them.
#
# Either way a¹ = (dξdx, dξdy, dξdz) and a² = (dηdx, dηdy, dηdz) are what the
# RHS wants:
#
#   ∇ₛf = a¹ ∂f/∂ξ + a² ∂f/∂η
#
# so a¹ = (dξdx, dξdy, dξdz) and a² = (dηdx, dηdy, dηdz), keeping Jexpresso's
# naming from metric_terms.jl with the third component added. The surface
# divergence of a Cartesian flux (F,G,H) is then
#
#   ∇ₛ·F = dξdx ∂F/∂ξ + dηdx ∂F/∂η
#        + dξdy ∂G/∂ξ + dηdy ∂G/∂η
#        + dξdz ∂H/∂ξ + dηdz ∂H/∂η
#
# which is the manifold version of the flat expansion in _expansion_inviscid!.
#
# Because a¹ and a² are tangential, the surface gradient of a scalar is
# automatically tangential — that is what makes the conservative pressure term
# of user_flux.jl produce no spurious normal force, as claimed there.
#
# THE MASS MATRIX
# ---------------
# Diagonal, by inexact (LGL) integration, assembled by direct stiffness
# summation over the elements that share a node:
#
#   M[ip] = Σ_e Σ_ij ω_i ω_j J[e,i,j]   over all (e,i,j) with connijk = ip
#
# Two properties make it a strong self-check, both reported by
# check_sphere_metrics:
#
#   * Σ_ip M[ip] is the SEM quadrature of the shell area and must equal 4πR²;
#   * aⁱ·a_j = δⁱⱼ must hold to round-off at every node.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export St_sphere_metrics
export build_sphere_metrics
export check_sphere_metrics
export sphere_metrics_form


struct St_sphere_metrics{TFloat}

    #  nelem × ngl × ngl
    Je   ::Array{TFloat, 3}     # surface Jacobian (the area element)

    dξdx ::Array{TFloat, 3}     # a¹ : the contravariant basis vector of ξ
    dξdy ::Array{TFloat, 3}
    dξdz ::Array{TFloat, 3}

    dηdx ::Array{TFloat, 3}     # a² : the contravariant basis vector of η
    dηdy ::Array{TFloat, 3}
    dηdz ::Array{TFloat, 3}

    nx   ::Array{TFloat, 3}     # outward UNIT normal n̂
    ny   ::Array{TFloat, 3}
    nz   ::Array{TFloat, 3}

    Jnx  ::Array{TFloat, 3}     # J n̂ : the AREA normal, i.e. the third metric
    Jny  ::Array{TFloat, 3}     # term J a³ of the radially extended map. It is
    Jnz  ::Array{TFloat, 3}     # what closes the metric identity (**), and it
                                # is NOT |J| times n̂ under :curl_invariant.

    M    ::Array{TFloat, 1}     # npoin : diagonal global mass matrix
    Minv ::Array{TFloat, 1}     # npoin : 1/M

    ξ    ::Array{TFloat, 1}     # ngl       LGL abscissae
    ω    ::Array{TFloat, 1}     # ngl       LGL weights
    dψ   ::Array{TFloat, 2}     # ngl × ngl dψ[k,i] = dψ_k/dξ at node i

    Δmin ::TFloat               # smallest distance between adjacent LGL nodes

    form ::Symbol               # :curl_invariant | :cross_product
end


"""
    sphere_metrics_form(inputs) -> Symbol

Which form of the manifold metric terms the case asked for, `:curl_invariant`
(the default, Kopriva 2006 Eq. (15) / Kelly et al. Eq. (23)) or
`:cross_product` (Kelly et al. Eq. (21)). Strings and the short names "ci"/"cp"
are accepted too, so that `:sphere_metrics => "CI"` in a user_inputs.jl does
what it looks like it does.
"""
function sphere_metrics_form(inputs)

    raw = get(inputs, :sphere_metrics, :curl_invariant)
    key = Symbol(lowercase(string(raw)))

    if key in (:curl_invariant, :curlinvariant, :ci, :curl, :invariant)
        return :curl_invariant
    elseif key in (:cross_product, :crossproduct, :cp, :cross)
        return :cross_product
    else
        error(string(" # ERROR sphere_metrics.jl: :sphere_metrics => ", repr(raw),
                     " is not a known form. Use :curl_invariant or :cross_product."))
    end
end


function build_sphere_metrics(mesh::St_mesh,
                              inputs;
                              backend = CPU(),
                              verbose = true,
                              TF      = TFloat)

    crd = mesh.coords          # (x,y,z); mesh.x/y/z are deprecated
    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    npoin = Int(mesh.npoin)
    form  = sphere_metrics_form(inputs)

    verbose && println(" # ")
    verbose && println(" # SPHERICAL SHELL METRICS .....................................")
    verbose && println(" #   metric terms: ",
                       form === :curl_invariant ?
                       "CURL-INVARIANT (Kopriva 2006 Eq. 15; Kelly et al. Eq. 23)" :
                       "CROSS-PRODUCT (Kelly et al. Eq. 21)")

    #-----------------------------------------------------------------------------
    # LGL points, weights and the differentiation matrix — the repo's own
    # routines, so the basis is bit-for-bit the one the flat cases use.
    # dψ[k,i] = dψ_k/dξ at node i (collocation: the quadrature points ARE the
    # interpolation points), which is the convention of _expansion_inviscid!.
    #-----------------------------------------------------------------------------
    lgl   = basis_structs_ξ_ω!(inputs[:interpolation_nodes], Int(mesh.nop), backend)
    ξ     = TF.(collect(lgl.ξ))
    ω     = TF.(collect(lgl.ω))
    basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξ, TF, backend)
    dψ    = TF.(Array(basis.dψ))

    Je   = zeros(TF, nelem, ngl, ngl)
    dξdx = zeros(TF, nelem, ngl, ngl); dξdy = zeros(TF, nelem, ngl, ngl); dξdz = zeros(TF, nelem, ngl, ngl)
    dηdx = zeros(TF, nelem, ngl, ngl); dηdy = zeros(TF, nelem, ngl, ngl); dηdz = zeros(TF, nelem, ngl, ngl)
    nx   = zeros(TF, nelem, ngl, ngl); ny   = zeros(TF, nelem, ngl, ngl); nz   = zeros(TF, nelem, ngl, ngl)
    Jnx  = zeros(TF, nelem, ngl, ngl); Jny  = zeros(TF, nelem, ngl, ngl); Jnz  = zeros(TF, nelem, ngl, ngl)
    M    = zeros(TF, npoin)

    # Typed function barrier: St_mesh fields carry no type annotation, so
    # mesh.connijk / mesh.coords read inside the node loops would box on every
    # access (the same reason the RHS kernels take their arrays as arguments).
    if form === :curl_invariant
        _sphere_metrics_curl_invariant!(Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                        nx, ny, nz, Jnx, Jny, Jnz, M,
                                        crd, mesh.connijk, dψ, ω,
                                        nelem, ngl, TF(mesh.radius))
    else
        _sphere_metrics_cross_product!(Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                       nx, ny, nz, Jnx, Jny, Jnz, M,
                                       crd, mesh.connijk, dψ, ω,
                                       nelem, ngl)
    end

    minimum(M) > 0 || error(" # ERROR sphere_metrics.jl: a node has zero mass — the grid is not covered by the elements.")
    Minv = one(TF) ./ M

    #-----------------------------------------------------------------------------
    # Smallest distance between adjacent LGL nodes, for the CFL condition. LGL
    # points cluster towards the element edges like 1/nop², so this is much
    # smaller than the element size and is what actually sets the time step.
    #-----------------------------------------------------------------------------
    Δmin = TF(Inf)
    for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl-1
            ip = mesh.connijk[iel,i,j]; iq = mesh.connijk[iel,i+1,j]
            Δmin = min(Δmin, sqrt((crd[1, ip]-crd[1, iq])^2 + (crd[2, ip]-crd[2, iq])^2 + (crd[3, ip]-crd[3, iq])^2))
        end
        for j = 1:ngl-1, i = 1:ngl
            ip = mesh.connijk[iel,i,j]; iq = mesh.connijk[iel,i,j+1]
            Δmin = min(Δmin, sqrt((crd[1, ip]-crd[1, iq])^2 + (crd[2, ip]-crd[2, iq])^2 + (crd[3, ip]-crd[3, iq])^2))
        end
    end

    metrics = St_sphere_metrics{TF}(Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                    nx, ny, nz, Jnx, Jny, Jnz,
                                    M, Minv, ξ, ω, dψ, Δmin, form)

    verbose && @printf(" #   smallest LGL node spacing Δmin = %.4e m\n", Δmin)
    verbose && println(" # SPHERICAL SHELL METRICS ..................................... DONE")

    return metrics
end


#---------------------------------------------------------------------------------
# _sphere_metrics_cross_product!  —  Kelly et al. Eq. (21).
#
#   n = a_ξ × a_η ,  J = |n| ,  n̂ = n/J ,  J n̂ = n
#   a¹ = (a_η × n)/J² ,  a² = (n × a_ξ)/J²
#
# aⁱ·a_j = δⁱⱼ and aⁱ·n̂ = 0 hold to round-off; the metric identity (**) holds
# only to the order of the polynomial approximation of the sphere.
#---------------------------------------------------------------------------------
function _sphere_metrics_cross_product!(Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                        nx, ny, nz, Jnx, Jny, Jnz,
                                        M::AbstractVector{TF},
                                        crd::AbstractMatrix{TF}, connijk,
                                        dψ::AbstractMatrix{TF}, ω::AbstractVector{TF},
                                        nelem::Int, ngl::Int) where {TF}

    @inbounds for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl

            a1x, a1y, a1z, a2x, a2y, a2z =
                _sphere_tangents(crd, connijk, dψ, iel, i, j, ngl, TF)

            # n = a_ξ × a_η ; J = |n|
            vx = a1y*a2z - a1z*a2y
            vy = a1z*a2x - a1x*a2z
            vz = a1x*a2y - a1y*a2x
            J  = sqrt(vx*vx + vy*vy + vz*vz)

            J > 0 || error(string(" # ERROR sphere_metrics.jl: zero surface Jacobian in element ", iel,
                                  " at node (", i, ",", j, "). The element is degenerate."))

            Je[iel,i,j]  = J
            nx[iel,i,j]  = vx/J; ny[iel,i,j]  = vy/J; nz[iel,i,j]  = vz/J
            Jnx[iel,i,j] = vx;   Jny[iel,i,j] = vy;   Jnz[iel,i,j] = vz

            invJ2 = one(TF)/(J*J)

            # a¹ = (a_η × n)/J²
            dξdx[iel,i,j] = (a2y*vz - a2z*vy)*invJ2
            dξdy[iel,i,j] = (a2z*vx - a2x*vz)*invJ2
            dξdz[iel,i,j] = (a2x*vy - a2y*vx)*invJ2

            # a² = (n × a_ξ)/J²
            dηdx[iel,i,j] = (vy*a1z - vz*a1y)*invJ2
            dηdy[iel,i,j] = (vz*a1x - vx*a1z)*invJ2
            dηdz[iel,i,j] = (vx*a1y - vy*a1x)*invJ2

            # direct stiffness summation of the diagonal mass matrix
            M[connijk[iel,i,j]] += ω[i]*ω[j]*J
        end
    end

    return nothing
end


#---------------------------------------------------------------------------------
# _sphere_metrics_curl_invariant!  —  Kopriva (2006) Eq. (15), Kelly et al.
# Eq. (23), specialised to the shell by the radial extension x → (ζ/R)x of the
# file header. Written out, with x̂ = x/|x| the unit radial at the node:
#
#   A_ξ = a_ξ × x̂ ,  A_η = a_η × x̂                (pointwise, nodal)
#
#   J a¹ =  A_η ,   J a² = -A_ξ ,   J = a_ξ·A_η = x̂·(a_ξ × a_η)
#
#   J n̂ = (R/2) [ ∂_η A_ξ - ∂_ξ A_η ]              (the surviving curl)
#
# The last line is the one that has to be DIFFERENTIATED, so it needs A_ξ and
# A_η at every node of the element before it can be formed: hence the two
# passes below. Because the same dψ differentiates both, the identity
#
#   ∂_ξ(J a¹) + ∂_η(J a²) + (2/R)(J n̂)
#     = ∂_ξ A_η - ∂_η A_ξ + [∂_η A_ξ - ∂_ξ A_η] = 0
#
# is an algebraic cancellation of the discrete derivatives — no smoothness of
# the grid, no accuracy of the interpolant, and no property of the sphere is
# used, which is exactly why it survives at finite order.
#
# Note also that J a¹ and J a² are cross products with x̂, so aⁱ·x̂ = 0 exactly
# (the basis is tangent to the TRUE sphere), and J is chosen as x̂·(a_ξ × a_η)
# rather than |a_ξ × a_η| so that aⁱ·a_j = δⁱⱼ is exact as well. n̂ is then x̂
# itself: on a shell the exact unit normal is known, and it is a better normal
# than any discrete cross product. The curl-form AREA normal J n̂ is a genuinely
# different vector and is kept separately, in Jnx/Jny/Jnz.
#---------------------------------------------------------------------------------
function _sphere_metrics_curl_invariant!(Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                         nx, ny, nz, Jnx, Jny, Jnz,
                                         M::AbstractVector{TF},
                                         crd::AbstractMatrix{TF}, connijk,
                                         dψ::AbstractMatrix{TF}, ω::AbstractVector{TF},
                                         nelem::Int, ngl::Int, R::TF) where {TF}

    Aξ = zeros(TF, ngl, ngl, 3)      # a_ξ × x̂ at every node of the element
    Aη = zeros(TF, ngl, ngl, 3)      # a_η × x̂

    halfR = TF(0.5)*R

    @inbounds for iel = 1:nelem

        #--- pass 1: the nodal products, and everything that follows from them
        for j = 1:ngl, i = 1:ngl

            a1x, a1y, a1z, a2x, a2y, a2z =
                _sphere_tangents(crd, connijk, dψ, iel, i, j, ngl, TF)

            ip = connijk[iel,i,j]
            xx, xy, xz = crd[1,ip], crd[2,ip], crd[3,ip]
            r  = sqrt(xx*xx + xy*xy + xz*xz)
            r > 0 || error(string(" # ERROR sphere_metrics.jl: node ", ip,
                                  " sits at the centre of the sphere."))
            xx /= r; xy /= r; xz /= r                       # x̂

            b1x = a1y*xz - a1z*xy                           # A_ξ = a_ξ × x̂
            b1y = a1z*xx - a1x*xz
            b1z = a1x*xy - a1y*xx

            b2x = a2y*xz - a2z*xy                           # A_η = a_η × x̂
            b2y = a2z*xx - a2x*xz
            b2z = a2x*xy - a2y*xx

            Aξ[i,j,1] = b1x; Aξ[i,j,2] = b1y; Aξ[i,j,3] = b1z
            Aη[i,j,1] = b2x; Aη[i,j,2] = b2y; Aη[i,j,3] = b2z

            # J = a_ξ·(a_η × x̂) = x̂·(a_ξ × a_η)
            J = a1x*b2x + a1y*b2y + a1z*b2z

            J > 0 || error(string(" # ERROR sphere_metrics.jl: non-positive surface Jacobian in element ",
                                  iel, " at node (", i, ",", j,
                                  "). The element is degenerate or wound inward."))

            Je[iel,i,j] = J
            nx[iel,i,j] = xx; ny[iel,i,j] = xy; nz[iel,i,j] = xz

            invJ = one(TF)/J
            dξdx[iel,i,j] =  b2x*invJ; dξdy[iel,i,j] =  b2y*invJ; dξdz[iel,i,j] =  b2z*invJ
            dηdx[iel,i,j] = -b1x*invJ; dηdy[iel,i,j] = -b1y*invJ; dηdz[iel,i,j] = -b1z*invJ

            M[ip] += ω[i]*ω[j]*J
        end

        #--- pass 2: the curl, J n̂ = (R/2)[∂_η A_ξ - ∂_ξ A_η]
        for j = 1:ngl, i = 1:ngl
            cx = cy = cz = zero(TF)
            for k = 1:ngl
                dη = dψ[k,j]
                dξ = dψ[k,i]
                cx += dη*Aξ[i,k,1] - dξ*Aη[k,j,1]
                cy += dη*Aξ[i,k,2] - dξ*Aη[k,j,2]
                cz += dη*Aξ[i,k,3] - dξ*Aη[k,j,3]
            end
            Jnx[iel,i,j] = halfR*cx
            Jny[iel,i,j] = halfR*cy
            Jnz[iel,i,j] = halfR*cz
        end
    end

    return nothing
end


#
# a_ξ = ∂x/∂ξ and a_η = ∂x/∂η at node (i,j) of element iel, from the nodal
# coordinates and the LGL differentiation matrix. dψ[k,i] = dψ_k/dξ at node i,
# which is the convention of _expansion_inviscid! and of the RHS kernels.
#
@inline function _sphere_tangents(crd::AbstractMatrix{TF}, connijk,
                                  dψ::AbstractMatrix{TF},
                                  iel::Int, i::Int, j::Int, ngl::Int, ::Type{TF}) where {TF}

    a1x = a1y = a1z = zero(TF)
    a2x = a2y = a2z = zero(TF)

    @inbounds for k = 1:ngl
        ipk = connijk[iel, k, j]
        dξ  = dψ[k,i]
        a1x += dξ*crd[1, ipk]; a1y += dξ*crd[2, ipk]; a1z += dξ*crd[3, ipk]

        ipl = connijk[iel, i, k]
        dη  = dψ[k,j]
        a2x += dη*crd[1, ipl]; a2y += dη*crd[2, ipl]; a2z += dη*crd[3, ipl]
    end

    return a1x, a1y, a1z, a2x, a2y, a2z
end


#---------------------------------------------------------------------------------
# check_sphere_metrics(mesh, metrics; …)
#
#   M1  metric identities .... aⁱ·a_j = δⁱⱼ at every node. This is the direct
#                              statement that the contravariant basis inverts
#                              the covariant one; it fails if the cross-product
#                              formulas or the differentiation matrix are wrong.
#   M2  tangency ............. a¹ and a² are perpendicular to n̂, so ∇ₛ of a
#                              scalar is tangential — the property the
#                              conservative pressure term of user_flux.jl relies
#                              on to produce no spurious normal force. Exact
#                              under both forms, but they mean slightly
#                              different things by it: :cross_product is
#                              tangent to the discrete surface (n̂ is its own
#                              cross product), :curl_invariant is tangent to the
#                              TRUE sphere (n̂ = x̂).
#   M3  outward normal ....... the AREA normal J n̂ points straight up,
#                              (J n̂/|J n̂|)·x̂ = 1 (the element orientation from
#                              orient_elements_outward! survived into the
#                              metrics). Under :cross_product J n̂ is J times the
#                              stored unit normal, so this is the same check it
#                              always was; under :curl_invariant the stored n̂ is
#                              x̂ by construction and it is the CURL-form area
#                              normal that has something to say here.
#   M4  mass matrix .......... Σ M[ip] is the SEM quadrature of the shell area
#                              and must equal 4πR². This one check exercises the
#                              Jacobian, the weights AND the direct stiffness
#                              summation at once.
#   M5  metric identity ...... ∂_ξ(J a¹) + ∂_η(J a²) + (2/R)(J n̂) = 0, the
#                              curved-surface form of the identity that makes a
#                              scheme constant-state preserving. Round-off under
#                              :curl_invariant — that is what the form is FOR —
#                              and O(hᴺ) under :cross_product.
#---------------------------------------------------------------------------------
function check_sphere_metrics(mesh::St_mesh, metrics::St_sphere_metrics;
                              verbose = true, atol_area = 1.0e-6, atol_normal = 1.0e-6,
                              atol_curvature = nothing)

    crd = mesh.coords          # (x,y,z); mesh.x/y/z are deprecated
    allok = true
    line(name, ok, extra="") = begin
        verbose && println("     ", ok ? "PASS" : "FAIL", "  ", name, extra == "" ? "" : string("  [", extra, "]"))
        ok
    end

    verbose && println(" # CHECK SPHERICAL SHELL METRICS ...............................")
    verbose && println("     metric terms: ", metrics.form)

    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    dψ    = metrics.dψ

    worst_id  = 0.0
    worst_tan = 0.0
    worst_nrm = 0.0

    for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl

            a1x = a1y = a1z = 0.0
            a2x = a2y = a2z = 0.0
            for k = 1:ngl
                ipk = mesh.connijk[iel, k, j]
                a1x += dψ[k,i]*crd[1, ipk]; a1y += dψ[k,i]*crd[2, ipk]; a1z += dψ[k,i]*crd[3, ipk]
                ipl = mesh.connijk[iel, i, k]
                a2x += dψ[k,j]*crd[1, ipl]; a2y += dψ[k,j]*crd[2, ipl]; a2z += dψ[k,j]*crd[3, ipl]
            end

            c1x, c1y, c1z = metrics.dξdx[iel,i,j], metrics.dξdy[iel,i,j], metrics.dξdz[iel,i,j]
            c2x, c2y, c2z = metrics.dηdx[iel,i,j], metrics.dηdy[iel,i,j], metrics.dηdz[iel,i,j]

            worst_id = max(worst_id,
                           abs(c1x*a1x + c1y*a1y + c1z*a1z - 1.0),
                           abs(c2x*a2x + c2y*a2y + c2z*a2z - 1.0),
                           abs(c1x*a2x + c1y*a2y + c1z*a2z),
                           abs(c2x*a1x + c2y*a1y + c2z*a1z))

            nxi, nyi, nzi = metrics.nx[iel,i,j], metrics.ny[iel,i,j], metrics.nz[iel,i,j]
            worst_tan = max(worst_tan,
                            abs(c1x*nxi + c1y*nyi + c1z*nzi)*mesh.radius,
                            abs(c2x*nxi + c2y*nyi + c2z*nzi)*mesh.radius)

            ip = mesh.connijk[iel,i,j]
            rr = sqrt(crd[1, ip]^2 + crd[2, ip]^2 + crd[3, ip]^2)

            # M3 reads the AREA normal J n̂, not the stored unit normal: under
            # :curl_invariant the latter IS x̂ and the check would be vacuous.
            Jnxi, Jnyi, Jnzi = metrics.Jnx[iel,i,j], metrics.Jny[iel,i,j], metrics.Jnz[iel,i,j]
            Jnmag = sqrt(Jnxi^2 + Jnyi^2 + Jnzi^2)
            worst_nrm = max(worst_nrm,
                            abs((Jnxi*crd[1, ip] + Jnyi*crd[2, ip] + Jnzi*crd[3, ip])/(Jnmag*rr) - 1.0))
        end
    end

    allok &= line("M1 metric identities aⁱ·a_j = δⁱⱼ", worst_id < 1.0e-10, @sprintf("max err = %.3e", worst_id))
    allok &= line("M2 contravariant basis is tangential", worst_tan < 1.0e-8, @sprintf("max |aⁱ·n̂|·R = %.3e", worst_tan))
    # M3 and M4 sit at round-off (1e-14) for nop >= 4 but are order-dependent
    # below that — the discrete tangents are derivatives of the polynomial
    # interpolant of the sphere, so the discrete normal and the discrete area
    # carry the interpolation error (measured 2.2e-8 and 2.7e-8 at nop=3 on a
    # 6x6-per-panel grid). The defaults are loose enough to hold at low order
    # and still tight enough that a wrong normal or a wrong mass matrix, both
    # O(1) errors, cannot pass.
    allok &= line("M3 outward area normal (J n̂)·x̂ = |J n̂|", worst_nrm < atol_normal,
                  @sprintf("max err = %.3e", worst_nrm))

    area  = sum(metrics.M)
    aref  = 4π*mesh.radius^2
    aerr  = abs(area - aref)/aref
    allok &= line("M4 Σ M[ip] = 4πR² (SEM quadrature of the area)", aerr < atol_area,
                  @sprintf("rel err = %.3e", aerr))

    #--- M5 the CURVED-SURFACE metric identity.
    #
    # On a FLAT grid the free-stream condition is ∇·(J aⁱ) = 0: the divergence
    # of a constant flux vanishes. On a CURVED surface it does NOT. Taking the
    # surface divergence of a constant vector c leaves the mean-curvature term
    #
    #   ∇ₛ·c = -H (c·n̂) ,   H = 2/R for a sphere of radius R,
    #
    # equivalently  ∂_ξ(J a¹) + ∂_η(J a²) = -(2/R)(J n̂), with J n̂ the AREA
    # normal held in Jnx/Jny/Jnz. That is the identity to test, and asserting
    # the flat one instead would flag a correct metric as broken (it reads
    # exactly 2/R).
    #
    # This is not a formality: it is the SAME curvature that produces
    # x·[∇ₛ·(φu⊗u)] = -φ|u|², which is what the Lagrange multiplier
    # μ = -φ|u|²/r² in user_source.jl cancels. If the discrete metrics got the
    # curvature wrong, the multiplier would not balance the discrete flux
    # divergence and the flow would drift off the shell.
    #
    # Evaluated with the same arithmetic the RHS uses, so it tests the real path.
    worst_fs = 0.0
    Hmean    = 2.0/mesh.radius
    for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl
            for comp = 1:3
                d = 0.0
                for k = 1:ngl
                    c1 = comp == 1 ? metrics.dξdx[iel,k,j] : comp == 2 ? metrics.dξdy[iel,k,j] : metrics.dξdz[iel,k,j]
                    c2 = comp == 1 ? metrics.dηdx[iel,i,k] : comp == 2 ? metrics.dηdy[iel,i,k] : metrics.dηdz[iel,i,k]
                    d += dψ[k,i]*metrics.Je[iel,k,j]*c1 + dψ[k,j]*metrics.Je[iel,i,k]*c2
                end
                Jnc = comp == 1 ? metrics.Jnx[iel,i,j] : comp == 2 ? metrics.Jny[iel,i,j] : metrics.Jnz[iel,i,j]
                worst_fs = max(worst_fs, abs(d + Hmean*Jnc)/(Hmean*metrics.Je[iel,i,j]))
            end
        end
    end
    #
    # THIS is where the two forms part company, and it is the reason
    # :curl_invariant is the default.
    #
    #   :curl_invariant — round-off, at every order and on any grid. J a¹, J a²
    #       and J n̂ are three parts of one discrete curl, so the identity is an
    #       algebraic cancellation between them rather than an approximation
    #       (Kopriva 2006, Sec. 7). It does not converge with nop because there
    #       is nothing to converge: measured on the shipped 10x10-per-panel
    #       cubed sphere, 1.6e-14 at nop=3, 3.9e-14 at nop=5, 6.8e-14 at nop=7,
    #       i.e. flat and at the level of the arithmetic. The tolerance below is
    #       1e-10: anything above that means the construction is broken, not
    #       that the grid is coarse.
    #
    #   :cross_product — resolution-dependent, because the discrete metric
    #       terms are polynomials of degree nop on a CURVED element, so the
    #       identity holds only to the order of the approximation. It converges
    #       spectrally with nop — same grid: 2.2e-2 at nop=3, 1.4e-4 at nop=5,
    #       5.4e-7 at nop=7 — which is twelve orders of magnitude worse than the
    #       curl form at nop=3 and still seven at nop=7. The tolerance is
    #       therefore loose on purpose: it is here to catch a WRONG identity — a
    #       sign slip or a missing term gives O(1), forty times the threshold —
    #       not to certify an accuracy.
    #
    # Pass atol_curvature explicitly to override either default.
    #
    atol_fs = atol_curvature === nothing ?
              (metrics.form === :curl_invariant ? 1.0e-10 : 5.0e-2) : atol_curvature
    allok &= line("M5 metric identity ∇ₛ·(J aⁱ) = -(2/R)(J n̂)", worst_fs < atol_fs,
                  @sprintf("rel err = %.3e, tol = %.1e (%s)", worst_fs, atol_fs,
                           metrics.form === :curl_invariant ? "exact by construction" :
                                                              "spectral in nop"))

    verbose && println(" # CHECK SPHERICAL SHELL METRICS ............................... ",
                       allok ? "ALL TESTS PASSED" : "FAILED")

    return allok
end


#---------------------------------------------------------------------------------
# Keep the momentum tangent to the shell.
#
# These two used to live in sphere_mesh.jl. They are not mesh construction —
# they act on the SOLUTION — so they moved here with the rest of the
# manifold-specific machinery when the grid builder was folded back into the
# ordinary gmsh reader in src/kernel/mesh/mesh.jl.
#
# `u` is the solution array indexed [ip, ivar], with the three Cartesian
# momentum (or velocity) components in columns ivar, ivar+1, ivar+2 — so
# ivar = 2 for the SWsphere state q = [φ, φu, φv, φw].
#---------------------------------------------------------------------------------

export project_momentum_to_sphere!
export sphere_normal_momentum

"""
    project_momentum_to_sphere!(u, mesh; ivar = 2) -> dmax

Remove the radial component of the momentum at every node, and return the
largest component removed, max|(φu)·x̂| — a drift diagnostic: watch it grow and
you are watching the discretization leave the manifold.
"""
function project_momentum_to_sphere!(u::AbstractMatrix, mesh::St_mesh; ivar::Int = 2)
    # Typed function barrier: St_mesh fields are ::Any (Base.@kwdef with no
    # annotations), so reading mesh.coords[ip,1] in the node loop would box on every
    # access. This runs four times per time step.
    return _project_momentum_kernel!(u, mesh.coords, Int(mesh.npoin), ivar)
end

function _project_momentum_kernel!(u::AbstractMatrix{TF}, crd::AbstractMatrix{TF},
                                   npoin::Int, ivar::Int) where {TF}

    dmax = zero(TF)

    @inbounds for ip = 1:npoin

        x, y, z = crd[1, ip], crd[2, ip], crd[3, ip]
        r2      = x*x + y*y + z*z

        mx = u[ip, ivar]
        my = u[ip, ivar+1]
        mz = u[ip, ivar+2]

        # μ = -(m·x)/r² ; the normal component removed is (m·x)/r
        mdotx = mx*x + my*y + mz*z
        μ     = -mdotx/r2

        u[ip, ivar]   = mx + μ*x
        u[ip, ivar+1] = my + μ*y
        u[ip, ivar+2] = mz + μ*z

        dmax = max(dmax, abs(mdotx)/sqrt(r2))
    end

    return dmax
end

"""
    sphere_normal_momentum(u, mesh; ivar = 2) -> dmax

The largest normal (off-shell) momentum component in the field, max|(φu)·x̂|,
WITHOUT modifying anything. `project_momentum_to_sphere!` returns the same
quantity as it removes it.
"""
function sphere_normal_momentum(u::AbstractMatrix, mesh::St_mesh; ivar::Int = 2)
    return _sphere_normal_momentum(u, mesh.coords, Int(mesh.npoin), ivar)
end

function _sphere_normal_momentum(u::AbstractMatrix{TF}, crd::AbstractMatrix{TF},
                                 npoin::Int, ivar::Int) where {TF}
    dmax = zero(TF)
    @inbounds for ip = 1:npoin
        x, y, z = crd[1, ip], crd[2, ip], crd[3, ip]
        dmax = max(dmax, abs(u[ip, ivar]*x + u[ip, ivar+1]*y + u[ip, ivar+2]*z) /
                         sqrt(x*x + y*y + z*z))
    end
    return dmax
end
