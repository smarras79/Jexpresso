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
# From those the contravariant basis is built. Kopriva, "Metric identities and
# the discontinuous spectral element method on curvilinear meshes", J. Sci.
# Comput. 26(3), 301-327 (2006), Eq. (15) — Sec. 3.2.3 / Eq. (23) of Kelly,
# Alves, Eckermann et al., JCP 552 (2026) 114683 — writes the metric terms of a
# 3-D map as a CURL,
#
#   J aⁱ = ½ [ ∂_ξk (∂x/∂ξʲ × x) − ∂_ξj (∂x/∂ξᵏ × x) ] ,  (i,j,k) cyclic,
#
# which is divergence-free by inspection, so the metric identities
# Σᵢ ∂_ξi (J aⁱ) = 0 survive DISCRETELY — the discrete divergence of a discrete
# curl vanishes identically when the same differentiation matrix is used for
# both. In 3-D that is a real constraint and it is what makes a scheme
# constant-state preserving.
#
# WHAT THAT FORM DEGENERATES TO ON A 2-D MANIFOLD  (read this before changing
# :sphere_metrics — it is not what one expects)
# ---------------------------------------------------------------------------
# A surface has no third reference coordinate to take the curl in, so one has
# to be invented: extend the surface off itself in some direction v(ξ,η),
#
#   X(ξ,η,ζ) = x(ξ,η) + (ζ − ζ₀) v(ξ,η) ,
#
# and evaluate Kopriva's formula on ζ = ζ₀. Doing that for i = 1 and i = 2, the
# ∂_η(v × x) terms cancel between the two halves of the curl and what is left is
#
#   J a¹ = a_η × v ,   J a² = v × a_ξ ,   J = v·(a_ξ × a_η)          (†)
#
# for ANY v. The curl structure imposes NOTHING on the two in-surface metric
# terms: (†) is the ordinary cross-product formula with v standing in for the
# surface normal. Choosing v = n̂ = (a_ξ × a_η)/|a_ξ × a_η| gives back exactly
# the textbook cross-product metrics — so on a 2-D manifold the cross-product
# form ALREADY IS the curl-invariant form, for the natural extension.
#
# The only object the curl does produce is the third metric term,
#
#   J a³ = J n̂ = (R/2) [ ∂_η (a_ξ × v) − ∂_ξ (a_η × v) ] ,          (‡)
#
# and by construction (†) and (‡) satisfy the curved-surface metric identity
#
#   ∂_ξ(J a¹) + ∂_η(J a²) + (2/R)(J a³) = 0                          (**)
#
# to round-off — again for ANY v (1.2e-16 measured directly on (†)/(‡), 1e-14
# once J a¹ is reassembled as J·a¹ the way the RHS stores it, at nop = 3, 5 and
# 7 under both choices below). The 2/R is the mean curvature of the sphere;
# (**), not the flat ∇·(J aⁱ) = 0, is the identity that holds on a curved
# surface.
#
# So (**) cannot be used to rank the two choices, and J a³ never enters the
# surface RHS. This is worth stating plainly because it is easy to get wrong:
# comparing (**) computed with J a³ from (‡) against (**) computed with
# J a³ = a_ξ × a_η measures the DEFINITION of the third term, not the quality
# of the metrics.
#
# WHAT ACTUALLY DIFFERS, AND THE TWO CHOICES OFFERED
# --------------------------------------------------
# What (†) leaves free is the extension direction v, and that choice is not
# empty — it decides which exact discrete properties hold.
#
#   :sphere_metrics => :cross_product   (default)   v = n̂, the discrete normal
#   :sphere_metrics => :radial                      v = x̂ = x/|x|, exact radial
#
#   v = n̂  · aⁱ·n̂ = 0 and J = |a_ξ × a_η| is the true area element of the
#           discrete surface;
#         · and, uniquely, the strong-form divergence ANNIHILATES A RIGID
#           ROTATION exactly. With u = Ω × x, D_ξu = Ω × a_ξ exactly, so
#             ∇ₛ·u = Ω·(a_ξ × a¹) + Ω·(a_η × a²)
#                  = Ω·[a_η(a_ξ·v) − a_ξ(a_η·v)]/J ,
#           which vanishes identically iff v ⊥ a_ξ, a_η, i.e. iff v ∥ a_ξ × a_η.
#           This is the manifold analogue of free-stream preservation for the
#           chain-rule form the RHS uses, and no other v has it.
#
#   v = x̂  · aⁱ·x̂ = 0, so ∇ₛ of a scalar is tangent to the TRUE sphere rather
#           than to the polynomial interpolant of it (radial leakage of the
#           surface gradient measured 2e-16 against 3e-7 for v = n̂ at nop=5) —
#           which is the property the Lagrange projection in user_source.jl is
#           there to enforce;
#         · but the rigid rotation is then only annihilated to O(hᴺ): for a
#           40 m/s rotation, 2.5e-9 / 6.0e-12 / 1.3e-14 s⁻¹ at nop = 3 / 5 / 7,
#           against 1e-19 flat for v = n̂.
#
# Measured on the shipped cubed sphere, the two are otherwise a wash: the
# surface divergence of the Galewsky jet and of φu agree to four significant
# figures, both dominated by the interpolation error of the jet itself
# (2.3e-2 relative at nop=5), and 6-day Galewsky runs finish with identical
# δE/E = -1.460e-03 and δmass/mass ~ 6e-12. The default is v = n̂ because it is
# the extension that keeps the exact property, and because it is what every
# result in this repository was produced with.
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

Which extension direction `v` of Eq. (†) in the header the case asked for:

  * `:cross_product` (default) — `v = n̂`, the discrete surface normal. This is
    Kopriva's curl-invariant form with the natural extension, and it reduces to
    the textbook cross-product metrics.
  * `:radial` — `v = x̂`, the exact radial direction.

`:curl_invariant` / `:ci` are accepted as aliases of `:radial`: on a 2-D
manifold BOTH choices are curl-invariant (see the header), so that name can
only mean "the one that is not the cross product". `:normal` / `:cp` are
aliases of `:cross_product`. Strings work too.
"""
function sphere_metrics_form(inputs)

    raw = get(inputs, :sphere_metrics, :cross_product)
    key = Symbol(lowercase(string(raw)))

    if key in (:cross_product, :crossproduct, :cp, :cross, :normal)
        return :cross_product
    elseif key in (:radial, :curl_invariant, :curlinvariant, :ci, :curl, :invariant)
        return :radial
    else
        error(string(" # ERROR sphere_metrics.jl: :sphere_metrics => ", repr(raw),
                     " is not a known form. Use :cross_product or :radial."))
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
    verbose && println(" #   Kopriva-form metric terms, extension direction v = ",
                       form === :radial ? "x̂ (radial)" : "n̂ (discrete normal; = the cross-product metrics)")

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
    _sphere_metrics!(Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                     nx, ny, nz, Jnx, Jny, Jnz, M,
                     crd, mesh.connijk, dψ, ω,
                     nelem, ngl, TF(mesh.radius), form === :radial)

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
# _sphere_metrics!  —  Eq. (†) and (‡) of the header, for either extension
# direction v:
#
#   A_ξ = a_ξ × v ,  A_η = a_η × v                      (pointwise, nodal)
#
#   J a¹ =  A_η ,   J a² = -A_ξ ,   J = a_ξ·A_η = v·(a_ξ × a_η)
#
#   J n̂ = (R/2) [ ∂_η A_ξ - ∂_ξ A_η ]                   (the curl third term)
#
# lradial = false  ⟹  v = n̂ = (a_ξ × a_η)/|a_ξ × a_η| — J is then |a_ξ × a_η|
#                     and J a¹ = (a_η × n)/J, i.e. the textbook cross-product
#                     metrics, unchanged from before :sphere_metrics existed.
# lradial = true   ⟹  v = x̂ = x/|x|.
#
# The last line is the one that has to be DIFFERENTIATED, so it needs A_ξ and
# A_η at every node of the element before it can be formed: hence the two
# passes. Because the same dψ differentiates both,
#
#   ∂_ξ(J a¹) + ∂_η(J a²) + (2/R)(J n̂)
#     = ∂_ξ A_η - ∂_η A_ξ + [∂_η A_ξ - ∂_ξ A_η] = 0
#
# identically — for EITHER v, which is the point made at length in the header:
# on a 2-D manifold that identity ranks nothing. What the choice of v does
# decide is which of the two exact properties holds, aⁱ·v = 0 with v the true
# normal (v = x̂) or the exact annihilation of a rigid rotation (v = n̂).
#---------------------------------------------------------------------------------
function _sphere_metrics!(Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                          nx, ny, nz, Jnx, Jny, Jnz,
                          M::AbstractVector{TF},
                          crd::AbstractMatrix{TF}, connijk,
                          dψ::AbstractMatrix{TF}, ω::AbstractVector{TF},
                          nelem::Int, ngl::Int, R::TF, lradial::Bool) where {TF}

    Aξ = zeros(TF, ngl, ngl, 3)      # a_ξ × v at every node of the element
    Aη = zeros(TF, ngl, ngl, 3)      # a_η × v

    halfR = TF(0.5)*R

    @inbounds for iel = 1:nelem

        #--- pass 1: the nodal products, and everything that follows from them
        for j = 1:ngl, i = 1:ngl

            a1x, a1y, a1z, a2x, a2y, a2z =
                _sphere_tangents(crd, connijk, dψ, iel, i, j, ngl, TF)

            ip = connijk[iel,i,j]

            if lradial
                vx, vy, vz = crd[1,ip], crd[2,ip], crd[3,ip]
            else
                vx = a1y*a2z - a1z*a2y
                vy = a1z*a2x - a1x*a2z
                vz = a1x*a2y - a1y*a2x
            end
            vn = sqrt(vx*vx + vy*vy + vz*vz)
            vn > 0 || error(string(" # ERROR sphere_metrics.jl: degenerate element ", iel,
                                   " at node (", i, ",", j,
                                   ") — zero surface normal, or a node at the centre of the sphere."))
            vx /= vn; vy /= vn; vz /= vn                     # v̂

            b1x = a1y*vz - a1z*vy                            # A_ξ = a_ξ × v
            b1y = a1z*vx - a1x*vz
            b1z = a1x*vy - a1y*vx

            b2x = a2y*vz - a2z*vy                            # A_η = a_η × v
            b2y = a2z*vx - a2x*vz
            b2z = a2x*vy - a2y*vx

            Aξ[i,j,1] = b1x; Aξ[i,j,2] = b1y; Aξ[i,j,3] = b1z
            Aη[i,j,1] = b2x; Aη[i,j,2] = b2y; Aη[i,j,3] = b2z

            # J = a_ξ·(a_η × v) = v·(a_ξ × a_η)
            J = a1x*b2x + a1y*b2y + a1z*b2z

            if !(J > 0)
                # A node that sits on a CUBE CORNER is the one degeneracy with a
                # known cause: the conformal cubed-sphere map has zero Jacobian
                # there and nothing can rescue it (THE 120° CORNER in
                # cubed_sphere_maps.jl). remap_cubed_sphere_nodes! refuses to
                # PRODUCE such a grid, but a grid READ with that map already
                # baked in arrives here, so name the suspect.
                x, y, z = crd[1,ip], crd[2,ip], crd[3,ip]
                r = sqrt(x*x + y*y + z*z)
                hi = max(abs(x), abs(y), abs(z)); lo = min(abs(x), abs(y), abs(z))
                oncorner = r > 0 && (hi - lo) < 1.0e-6*r
                error(string(" # ERROR sphere_metrics.jl: non-positive surface Jacobian in element ",
                             iel, " at node (", i, ",", j,
                             "). The element is degenerate or wound inward.",
                             oncorner ?
                             string("\n #   That node sits on a CUBE CORNER. If this grid carries the conformal",
                                    "\n #   cubed-sphere map, that is the map's own singularity — its Jacobian is",
                                    "\n #   zero at all eight corners — and it cannot be regularised away; see THE",
                                    "\n #   120° CORNER in src/kernel/mesh/cubed_sphere_maps.jl. Rebuild the grid",
                                    "\n #   with the equiangular map (tools/generate_cubed_sphere.jl).") : ""))
            end

            Je[iel,i,j] = J
            nx[iel,i,j] = vx; ny[iel,i,j] = vy; nz[iel,i,j] = vz

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
#                              metrics).
#   M4  mass matrix .......... Σ M[ip] is the SEM quadrature of the shell area
#                              and must equal 4πR². This one check exercises the
#                              Jacobian, the weights AND the direct stiffness
#                              summation at once.
#   M5  metric identity ...... ∂_ξ(J a¹) + ∂_η(J a²) + (2/R)(J n̂) = 0. Round-off
#                              under BOTH forms — it is an algebraic property of
#                              the curl definition of J n̂, so it checks THIS
#                              FILE, not the geometry, and it ranks nothing.
#   M6  geometry ............. the curl third term J a³ against the plain area
#                              normal a_ξ × a_η. This is the truncation error of
#                              the element geometry — the number M5 would report
#                              if J a³ were not defined by the curl.
#   M7  rigid rotation ....... ∇ₛ·(Ω × x) = 0 in the strong form. The manifold
#                              analogue of free-stream preservation, and the one
#                              check that really does separate the two forms:
#                              exact for :cross_product, O(hᴺ) for :radial.
#---------------------------------------------------------------------------------
function check_sphere_metrics(mesh::St_mesh, metrics::St_sphere_metrics;
                              verbose = true, atol_area = 1.0e-6, atol_normal = 1.0e-6,
                              atol_curvature = nothing, atol_geometry = 5.0e-2)

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
    worst_geo = 0.0

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
            # :radial the latter IS x̂ and the check would be vacuous.
            Jnxi, Jnyi, Jnzi = metrics.Jnx[iel,i,j], metrics.Jny[iel,i,j], metrics.Jnz[iel,i,j]
            Jnmag = sqrt(Jnxi^2 + Jnyi^2 + Jnzi^2)
            worst_nrm = max(worst_nrm,
                            abs((Jnxi*crd[1, ip] + Jnyi*crd[2, ip] + Jnzi*crd[3, ip])/(Jnmag*rr) - 1.0))

            # M6: the curl third term (‡) against the plain area normal
            # a_ξ × a_η. They agree only in the continuum, so this is the honest
            # measure of how well the element resolves the curvature — the
            # number that M5 would report if J a³ were NOT defined by the curl.
            worst_geo = max(worst_geo,
                            sqrt((Jnxi - (a1y*a2z - a1z*a2y))^2 +
                                 (Jnyi - (a1z*a2x - a1x*a2z))^2 +
                                 (Jnzi - (a1x*a2y - a1y*a2x))^2)/metrics.Je[iel,i,j])
        end
    end

    #--- M7: does the strong-form surface divergence annihilate a RIGID
    # ROTATION? u = Ω × x is tangential and divergence free on any sphere, and
    # D_ξu = Ω × a_ξ holds exactly for the nodal interpolant, so whatever comes
    # out is pure metric error. This is the manifold analogue of free-stream
    # preservation for the chain-rule form the RHS uses — the flat "constant
    # flux stays constant" test is empty here, because a constant flux has zero
    # discrete derivative whatever the metrics are.
    #
    # It is EXACT for :cross_product and only O(hᴺ) for :radial; see the header
    # for why v ∥ a_ξ × a_η is the unique choice that has it.
    worst_rot = 0.0
    Ω = (0.3, -0.7, 0.5) ./ sqrt(0.83) .* (40.0/mesh.radius)   # a 40 m/s rotation
    for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl
            d = 0.0
            for c = 1:3
                dξ = dη = 0.0
                for k = 1:ngl
                    ipk = mesh.connijk[iel,k,j]; ipl = mesh.connijk[iel,i,k]
                    uk = Ω[mod1(c+1,3)]*crd[mod1(c+2,3), ipk] - Ω[mod1(c+2,3)]*crd[mod1(c+1,3), ipk]
                    ul = Ω[mod1(c+1,3)]*crd[mod1(c+2,3), ipl] - Ω[mod1(c+2,3)]*crd[mod1(c+1,3), ipl]
                    dξ += dψ[k,i]*uk
                    dη += dψ[k,j]*ul
                end
                c1 = c == 1 ? metrics.dξdx[iel,i,j] : c == 2 ? metrics.dξdy[iel,i,j] : metrics.dξdz[iel,i,j]
                c2 = c == 1 ? metrics.dηdx[iel,i,j] : c == 2 ? metrics.dηdy[iel,i,j] : metrics.dηdz[iel,i,j]
                d += dξ*c1 + dη*c2
            end
            worst_rot = max(worst_rot, abs(d))
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
    #
    # THESE TOLERANCES ARE NOT MIS-CALIBRATED — the question comes up whenever a
    # bad grid fails them. On the shipped 10-per-panel equiangular cubed sphere
    # they hold with four to nine orders of margin at every order, :radial:
    #
    #   nop      3        4        5        6        7        8
    #   M3    1.5e-07  1.5e-10  1.3e-12  1.3e-15  5.6e-16  4.4e-16   (tol 1e-6)
    #   M4    2.2e-08  1.3e-11  4.4e-14  4.2e-14  2.7e-14  3.6e-14   (tol 1e-6)
    #   M6    4.5e-04  2.7e-05  1.2e-06  7.8e-08  2.9e-09  1.9e-10   (tol 5e-2)
    #   M7    2.5e-09  2.3e-12  6.0e-12  5.3e-15  1.3e-14  1.3e-17   (tol 1e-8)
    #
    # M6 falling by ~1.5 orders per order of nop is the spectral convergence a
    # SMOOTH map has to show. A map that fails these checks and does NOT converge
    # in nop is not being judged too harshly; it is not smooth. That is exactly
    # what happened to the corner-stretched conformal map (M6 = 0.41 flat in both
    # nop and h) — see THE 120° CORNER in cubed_sphere_maps.jl.
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
    # Round-off for BOTH forms, at every order and on any grid: J a¹, J a² and
    # J n̂ are three parts of one discrete curl, so this is an algebraic
    # cancellation between them and not an approximation (Kopriva 2006, Sec. 7,
    # specialised in the header). Measured 1.5e-14 / 3.9e-14 / 7.0e-14 at
    # nop = 3 / 5 / 7 under :cross_product and 1.6e-14 / 3.9e-14 / 6.8e-14 under
    # :radial — flat in nop, and at the level of the arithmetic (it is 1.2e-16
    # on the raw J aⁱ; the extra two digits are the J·(J aⁱ/J) round trip this
    # check has to make, since the metrics are stored split).
    #
    # It therefore CANNOT be used to rank the two forms — it is a check on this
    # file, not on the geometry. What it does catch, and catch hard, is a sign
    # slip or a missing term in (†)/(‡), which give O(1). The number that DOES
    # measure how well the element resolves the curvature is M6.
    #
    atol_fs = atol_curvature === nothing ? 1.0e-10 : atol_curvature
    allok &= line("M5 metric identity ∇ₛ·(J aⁱ) = -(2/R)(J n̂)", worst_fs < atol_fs,
                  @sprintf("rel err = %.3e, tol = %.1e (exact by construction)", worst_fs, atol_fs))

    #--- M6 the curl third term against the plain area normal a_ξ × a_η. These
    # agree only in the continuum, so this is the truncation error of the
    # element geometry, and it converges spectrally. Measured on the shipped
    # cubed sphere: :cross_product 2.2e-2 / 1.4e-4 / 5.4e-7 at nop = 3/5/7, and
    # :radial 3.4e-4 / 1.1e-6 / 2.8e-9. Loose on purpose — it is here to catch a
    # WRONG third term, which is O(1), not to certify an accuracy.
    allok &= line("M6 curl normal J a³ vs a_ξ × a_η", worst_geo < atol_geometry,
                  @sprintf("rel err = %.3e (spectral in nop)", worst_geo))

    #--- M7 rigid-rotation annihilation, computed above. Normalised by |Ω|,
    # i.e. by the reciprocal time scale of the rotation, so the number is the
    # spurious divergence a 40 m/s solid-body flow would feel, in s⁻¹.
    atol_rot = metrics.form === :cross_product ? 1.0e-16 : 1.0e-8
    allok &= line("M7 ∇ₛ·(Ω × x) = 0 (rigid rotation, 40 m/s)", worst_rot < atol_rot,
                  @sprintf("max = %.3e s⁻¹, tol = %.1e (%s)", worst_rot, atol_rot,
                           metrics.form === :cross_product ? "exact for v = n̂" :
                                                             "O(hᴺ) for v = x̂"))

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
