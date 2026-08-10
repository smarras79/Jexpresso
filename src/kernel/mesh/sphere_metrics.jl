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
# Their cross product is normal to the surface and its length is the surface
# Jacobian — the area element:
#
#   n = a_ξ × a_η ,   J = |n| ,   n̂ = n/J
#
# The contravariant basis is then
#
#   a¹ = (a_η × n)/J² ,   a² = (n × a_ξ)/J²
#
# which satisfies aⁱ·a_j = δⁱⱼ and, being cross products involving n, is
# TANGENT to the shell. Those are exactly the quantities the RHS wants:
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


struct St_sphere_metrics{TFloat}

    #  nelem × ngl × ngl
    Je   ::Array{TFloat, 3}     # surface Jacobian |a_ξ × a_η| (the area element)

    dξdx ::Array{TFloat, 3}     # a¹ : the contravariant basis vector of ξ
    dξdy ::Array{TFloat, 3}
    dξdz ::Array{TFloat, 3}

    dηdx ::Array{TFloat, 3}     # a² : the contravariant basis vector of η
    dηdy ::Array{TFloat, 3}
    dηdz ::Array{TFloat, 3}

    nx   ::Array{TFloat, 3}     # outward unit normal n̂ = (a_ξ × a_η)/J
    ny   ::Array{TFloat, 3}
    nz   ::Array{TFloat, 3}

    M    ::Array{TFloat, 1}     # npoin : diagonal global mass matrix
    Minv ::Array{TFloat, 1}     # npoin : 1/M

    ξ    ::Array{TFloat, 1}     # ngl       LGL abscissae
    ω    ::Array{TFloat, 1}     # ngl       LGL weights
    dψ   ::Array{TFloat, 2}     # ngl × ngl dψ[k,i] = dψ_k/dξ at node i

    Δmin ::TFloat               # smallest distance between adjacent LGL nodes
end


function build_sphere_metrics(mesh::St_mesh_sphere,
                              inputs;
                              backend = CPU(),
                              verbose = true,
                              TF      = TFloat)

    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    npoin = Int(mesh.npoin)

    verbose && println(" # ")
    verbose && println(" # SPHERICAL SHELL METRICS .....................................")

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
    M    = zeros(TF, npoin)

    for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl

            # covariant (tangent) vectors from the nodal coordinates
            a1x = a1y = a1z = zero(TF)
            a2x = a2y = a2z = zero(TF)
            @inbounds for k = 1:ngl
                ipk = mesh.connijk[iel, k, j]
                a1x += dψ[k,i]*mesh.x[ipk]; a1y += dψ[k,i]*mesh.y[ipk]; a1z += dψ[k,i]*mesh.z[ipk]

                ipl = mesh.connijk[iel, i, k]
                a2x += dψ[k,j]*mesh.x[ipl]; a2y += dψ[k,j]*mesh.y[ipl]; a2z += dψ[k,j]*mesh.z[ipl]
            end

            # n = a_ξ × a_η ; J = |n|
            vx = a1y*a2z - a1z*a2y
            vy = a1z*a2x - a1x*a2z
            vz = a1x*a2y - a1y*a2x
            J  = sqrt(vx*vx + vy*vy + vz*vz)

            J > 0 || error(string(" # ERROR sphere_metrics.jl: zero surface Jacobian in element ", iel,
                                  " at node (", i, ",", j, "). The element is degenerate."))

            Je[iel,i,j] = J
            nx[iel,i,j] = vx/J; ny[iel,i,j] = vy/J; nz[iel,i,j] = vz/J

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
            M[mesh.connijk[iel,i,j]] += ω[i]*ω[j]*J
        end
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
            Δmin = min(Δmin, sqrt((mesh.x[ip]-mesh.x[iq])^2 + (mesh.y[ip]-mesh.y[iq])^2 + (mesh.z[ip]-mesh.z[iq])^2))
        end
        for j = 1:ngl-1, i = 1:ngl
            ip = mesh.connijk[iel,i,j]; iq = mesh.connijk[iel,i,j+1]
            Δmin = min(Δmin, sqrt((mesh.x[ip]-mesh.x[iq])^2 + (mesh.y[ip]-mesh.y[iq])^2 + (mesh.z[ip]-mesh.z[iq])^2))
        end
    end

    metrics = St_sphere_metrics{TF}(Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                    nx, ny, nz, M, Minv, ξ, ω, dψ, Δmin)

    verbose && @printf(" #   smallest LGL node spacing Δmin = %.4e m\n", Δmin)
    verbose && println(" # SPHERICAL SHELL METRICS ..................................... DONE")

    return metrics
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
#                              on to produce no spurious normal force.
#   M3  outward normal ....... n̂·x̂ = 1 (the element orientation from
#                              orient_elements_outward! survived into the metrics)
#   M4  mass matrix .......... Σ M[ip] is the SEM quadrature of the shell area
#                              and must equal 4πR². This one check exercises the
#                              Jacobian, the weights AND the direct stiffness
#                              summation at once.
#   M5  constant-state ....... the surface divergence of a CONSTANT flux must
#                              vanish (the discrete "metric identity" / free-
#                              stream preservation). A scheme that fails this
#                              generates spurious motion out of a state at rest.
#---------------------------------------------------------------------------------
function check_sphere_metrics(mesh::St_mesh_sphere, metrics::St_sphere_metrics;
                              verbose = true, atol_area = 1.0e-6, atol_normal = 1.0e-6,
                              atol_curvature = 5.0e-2)

    allok = true
    line(name, ok, extra="") = begin
        verbose && println("     ", ok ? "PASS" : "FAIL", "  ", name, extra == "" ? "" : string("  [", extra, "]"))
        ok
    end

    verbose && println(" # CHECK SPHERICAL SHELL METRICS ...............................")

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
                a1x += dψ[k,i]*mesh.x[ipk]; a1y += dψ[k,i]*mesh.y[ipk]; a1z += dψ[k,i]*mesh.z[ipk]
                ipl = mesh.connijk[iel, i, k]
                a2x += dψ[k,j]*mesh.x[ipl]; a2y += dψ[k,j]*mesh.y[ipl]; a2z += dψ[k,j]*mesh.z[ipl]
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
            rr = sqrt(mesh.x[ip]^2 + mesh.y[ip]^2 + mesh.z[ip]^2)
            worst_nrm = max(worst_nrm,
                            abs((nxi*mesh.x[ip] + nyi*mesh.y[ip] + nzi*mesh.z[ip])/rr - 1.0))
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
    allok &= line("M3 outward unit normal n̂·x̂ = 1", worst_nrm < atol_normal, @sprintf("max err = %.3e", worst_nrm))

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
    # equivalently  ∂_ξ(J a¹) + ∂_η(J a²) = -(2J/R) n̂. That is the identity to
    # test, and asserting the flat one instead would flag a correct metric as
    # broken (it reads exactly 2/R).
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
                nc = comp == 1 ? metrics.nx[iel,i,j] : comp == 2 ? metrics.ny[iel,i,j] : metrics.nz[iel,i,j]
                worst_fs = max(worst_fs, abs(d/metrics.Je[iel,i,j] + Hmean*nc)/Hmean)
            end
        end
    end
    # Resolution-dependent, unlike M1-M4: the discrete metric terms are
    # polynomials of degree nop on a CURVED element, so this identity holds only
    # to the order of the approximation — and converges spectrally with nop
    # (measured on a 6x6-per-panel cubed sphere: 1.2e-5 at nop=4, 1.3e-8 at
    # nop=6, 2.4e-11 at nop=8; on a very coarse 3x3 panel it is 2.6e-2 at nop=3
    # and 3.3e-4 at nop=5). The default tolerance is therefore loose on purpose:
    # it is here to catch a WRONG identity — a sign slip or a missing term gives
    # O(1), forty times the threshold — not to certify an accuracy. Tighten it
    # with atol_curvature when you know the resolution.
    allok &= line("M5 curvature identity ∇ₛ·(J aⁱ)/J = -(2/R) n̂", worst_fs < atol_curvature,
                  @sprintf("rel err = %.3e (spectral in nop)", worst_fs))

    verbose && println(" # CHECK SPHERICAL SHELL METRICS ............................... ",
                       allok ? "ALL TESTS PASSED" : "FAILED")

    return allok
end
