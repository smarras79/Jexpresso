# =============================================================================
#  Element Learning — NON-CONSTANT DIFFUSIVITY  (Option 1: synthesized Jacobians)
# =============================================================================
#
#  Governing local problem on element K, written on the reference element
#  [-1,1]^d  (d = 2 here):
#
#     (A_{v^o,v})_{ij} = ∫_K a(x) ∇_x φ_i · ∇_x φ_j
#                      = ∫_{[-1,1]^2} (â ∇_ξ φ_i) · (∇_ξ φ_j) ,
#
#  where the *reference diffusivity* â is the 2×2 SPD tensor that fuses the
#  physical (scalar) diffusivity a and the element Jacobian J_K (= ∂x/∂ξ):
#
#     â(ξ) = (|K|/|[-1,1]^d|) · a · J_K^{-1} J_K^{-T}
#          = det(J_K) · a · (J_K^T J_K)^{-1}            (affine element).
#
#  The element-learning surrogate maps â ↦ T^{ie} = (A_{v^o,v^o})^{-1} A_{v^o,v^b}.
#
#  Feature representation (per the design decision): the 3 unique entries
#  (â11, â12, â22) of the SPD tensor â at each of the (k+1)^2 reference nodes,
#  i.e. a 3·(k+1)^2 input vector. For an AFFINE element (parallelogram) â is
#  ξ-independent, so the per-node entries are identical; for a general (bilinear)
#  straight-sided quad or a curved element â varies within the element, and the
#  per-node layout carries that variation.
#
#  Sampling source (Option 1): random element shapes are produced by
#  *synthesizing* random affine Jacobians J_K directly (no mesh file needed),
#  which induces the distribution of â matrices used to train the surrogate.
#
#  This whole path is gated by inputs[:lEL_nonconstant] (default false), so the
#  existing constant-amplitude pipeline and its trained model are untouched.
#
#  Flags (user_inputs.jl):
#    :lEL_nonconstant => true   activate the non-constant â feature in BOTH the
#                               sampler and inference (needs a model trained on
#                               the 3·(k+1)² feature, e.g. via tools/EL_training).
#    :lEL_xidependent => true   (sampling only) synthesize ξ-dependent â for
#                               curved / non-affine elements; default false uses
#                               affine Jacobians (constant â per element).
#
#  Sampling and inference share ONE node ordering (mesh.conn order, via
#  el_conn_to_ij), so the per-element â feature built from metrics at inference
#  matches the synthesized â feature used for training, and the recovery in
#  elementLearning_infer! is reused unchanged. The physical scalar diffusivity at
#  inference is el_diffusivity(x,y) (default 1; redefine for a varying coeff.).
# =============================================================================

using LinearAlgebra
using Random

# ── Reference-element node partition ─────────────────────────────────────────
# 2D tensor nodes (i,j), i,j ∈ 1:ngl, flattened P = (i-1)*ngl + j.
# Interior (v^o):  2 ≤ i ≤ ngl-1 AND 2 ≤ j ≤ ngl-1.
# Boundary (v^b):  everything on the reference perimeter (i or j ∈ {1,ngl}).
# Returns (vo, vb) as Vectors of flattened node indices in increasing P order.
function el_reference_node_partition(ngl::Int)
    vo = Int[]; vb = Int[]
    for i = 1:ngl, j = 1:ngl
        P = (i-1)*ngl + j
        if 2 <= i <= ngl-1 && 2 <= j <= ngl-1
            push!(vo, P)
        else
            push!(vb, P)
        end
    end
    return vo, vb
end

# ── conn ↔ tensor-(i,j) mapping ──────────────────────────────────────────────
# The element-learning recovery (elementLearning_infer!) addresses local nodes
# in mesh.conn order: boundary nodes occupy conn columns 1:elnbdypoints and
# element-interior nodes occupy elnbdypoints+1:nelpoints (conn is built
# corners→edges→volume, so the interior is last). The metrics and connijk are
# in tensor (i,j) order. To make the SAMPLED training data and the INFERENCE
# feature/recovery use one common ordering, we map every conn-local index m to
# its tensor position (i,j) once (the layout is identical for every element).
function el_conn_to_ij(mesh, ngl::Int)
    iel = 1
    gid2ij = Dict{Int,Tuple{Int,Int}}()
    for i = 1:ngl, j = 1:ngl
        gid2ij[mesh.connijk[iel, i, j, 1]] = (i, j)
    end
    npts    = ngl*ngl
    conn2ij = Vector{Tuple{Int,Int}}(undef, npts)
    for m = 1:npts
        conn2ij[m] = gid2ij[mesh.conn[iel, m]]
    end
    return conn2ij
end

# Flattened tensor index P = (i-1)*ngl + j for each conn-local node m.
el_conn_perm(conn2ij, ngl::Int) = [ (ij[1]-1)*ngl + ij[2] for ij in conn2ij ]

# ── Precompute the three reference-stiffness component matrices ───────────────
# Each is (ngl² × ngl²):
#   K11_{PQ} = Σ_{k,l} ω_k ω_l (∂_ξ φ_P)(∂_ξ φ_Q)
#   K12_{PQ} = Σ_{k,l} ω_k ω_l ((∂_ξ φ_P)(∂_η φ_Q) + (∂_η φ_P)(∂_ξ φ_Q))
#   K22_{PQ} = Σ_{k,l} ω_k ω_l (∂_η φ_P)(∂_η φ_Q)
# so that for a CONSTANT â the local stiffness is the linear combination
#   A_loc = â11·K11 + â12·K12 + â22·K22.
# Convention (matching build_laplace_matrix, NSD_2D): ψ[i,k]=ℓ_i(ξ_k),
# dψ[i,k]=ℓ'_i(ξ_k); 2D node P=(i,j) has ∂_ξφ_P|_{k,l}=dψ[i,k]ψ[j,l],
# ∂_ηφ_P|_{k,l}=ψ[i,k]dψ[j,l].
function el_precompute_reference_components(ψ::AbstractMatrix, dψ::AbstractMatrix,
                                            ω::AbstractVector, ngl::Int)
    n  = ngl*ngl
    K11 = zeros(n, n); K12 = zeros(n, n); K22 = zeros(n, n)
    for l = 1:ngl, k = 1:ngl
        w = ω[k]*ω[l]
        for i = 1:ngl, j = 1:ngl
            P   = (i-1)*ngl + j
            gξP = dψ[i,k]*ψ[j,l]
            gηP = ψ[i,k]*dψ[j,l]
            for m = 1:ngl, n2 = 1:ngl
                Q   = (m-1)*ngl + n2
                gξQ = dψ[m,k]*ψ[n2,l]
                gηQ = ψ[m,k]*dψ[n2,l]
                K11[P,Q] += w*gξP*gξQ
                K12[P,Q] += w*(gξP*gηQ + gηP*gξQ)
                K22[P,Q] += w*gηP*gηQ
            end
        end
    end
    return K11, K12, K22
end

# ── Assemble the local reference stiffness for a CONSTANT â ───────────────────
el_assemble_local_stiffness(K11, K12, K22, â11, â12, â22) =
    â11 .* K11 .+ â12 .* K12 .+ â22 .* K22

# ── Assemble the local reference stiffness for a ξ-DEPENDENT (per-node) â ─────
# `ahat_nodes` is (3, ngl²) in tensor flat order P=(i-1)*ngl+j (rows â11,â12,â22).
# Full quadrature assembly — used for curved / non-affine elements where â varies
# across the reference element. Returns the (ngl² × ngl²) stiffness in tensor
# node order.
function el_assemble_local_stiffness_field(ψ::AbstractMatrix, dψ::AbstractMatrix,
                                           ω::AbstractVector,
                                           ahat_nodes::AbstractMatrix, ngl::Int)
    n = ngl*ngl
    A = zeros(n, n)
    for l = 1:ngl, k = 1:ngl
        q   = (k-1)*ngl + l
        a11 = ahat_nodes[1,q]; a12 = ahat_nodes[2,q]; a22 = ahat_nodes[3,q]
        w   = ω[k]*ω[l]
        for i = 1:ngl, j = 1:ngl
            P   = (i-1)*ngl + j
            gξP = dψ[i,k]*ψ[j,l]
            gηP = ψ[i,k]*dψ[j,l]
            for m = 1:ngl, n2 = 1:ngl
                Q   = (m-1)*ngl + n2
                gξQ = dψ[m,k]*ψ[n2,l]
                gηQ = ψ[m,k]*dψ[n2,l]
                A[P,Q] += w*( gξP*(a11*gξQ + a12*gηQ) + gηP*(a12*gξQ + a22*gηQ) )
            end
        end
    end
    return A
end

# ── Synthesize a random affine Jacobian → SPD reference diffusivity â ─────────
# Draw an element shape J_K and a scalar amplitude a ~ U(amin,amax); return the
# constant SPD tensor â = a·det(J)·(JᵀJ)^{-1} as its 3 unique entries, plus J.
# The shape model: rotation θ, anisotropic stretch (s1,s2), and a shear γ — a
# generic distribution of well-formed (det>0) elements.
function synthesize_random_ahat(rng::AbstractRNG;
                                amin=1.0, amax=1.0,   # a is a nuisance (Tie is amplitude-invariant); a=1 matches inference
                                smin=0.5, smax=2.0, shear_max=0.5)
    a  = amin + (amax-amin)*rand(rng)
    θ  = 2π*rand(rng)
    s1 = smin + (smax-smin)*rand(rng)
    s2 = smin + (smax-smin)*rand(rng)
    γ  = shear_max*(2*rand(rng) - 1)
    R  = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    S  = [s1 γ; 0.0 s2]                      # stretch + shear (upper triangular)
    J  = R*S                                  # affine Jacobian J_K = ∂x/∂ξ
    â  = a * det(J) * inv(transpose(J)*J)     # = (|K|/|ref|)·a·J^{-1}J^{-T}, SPD
    return â[1,1], â[1,2], â[2,2], J
end

# ── ξ-DEPENDENT â: synthesize a per-node SPD field (curved / non-affine) ─────
# Draws a base affine Jacobian J0 and a smooth, anisotropic per-node warp so the
# induced node Jacobian J(ξ) — and hence â(ξ) = a·det(J)·(JᵀJ)⁻¹ — varies across
# the reference element. Normalized index coordinates ξ̂_i = -1+2(i-1)/(ngl-1) in
# [-1,1] drive a low-frequency warp (no external node coordinates needed). The
# warp magnitude `curve` is kept small and det(J) is guarded > 0. Returns a
# (3, ngl²) matrix in tensor flat order P=(i-1)*ngl+j: rows (â11, â12, â22).
function synthesize_random_ahat_field(rng::AbstractRNG, ngl::Int;
                                      amin=1.0, amax=1.0,   # a is a nuisance (Tie is amplitude-invariant); a=1 matches inference
                                      smin=0.5, smax=2.0, shear_max=0.5,
                                      curve=0.15)
    a  = amin + (amax-amin)*rand(rng)
    θ  = 2π*rand(rng)
    s1 = smin + (smax-smin)*rand(rng)
    s2 = smin + (smax-smin)*rand(rng)
    γ  = shear_max*(2*rand(rng) - 1)
    J0 = [cos(θ) -sin(θ); sin(θ) cos(θ)] * [s1 γ; 0.0 s2]

    # smooth anisotropic warp coefficients (small ⇒ det stays positive)
    c11 = curve*(2*rand(rng)-1); c12 = curve*(2*rand(rng)-1)
    c21 = curve*(2*rand(rng)-1); c22 = curve*(2*rand(rng)-1)

    coord(i) = -1.0 + 2.0*(i-1)/(ngl-1)
    ah = Matrix{Float64}(undef, 3, ngl*ngl)
    @inbounds for i = 1:ngl, j = 1:ngl
        P  = (i-1)*ngl + j
        ξ̂  = coord(i);  η̂ = coord(j)
        W  = [1.0 + c11*ξ̂   c12*η̂ ;  c21*ξ̂   1.0 + c22*η̂]   # per-node anisotropic warp
        J  = J0 * W
        d  = det(J)
        if d <= 0                                            # guard: keep det>0
            J = J0;  d = det(J0)
        end
        â  = a * d * inv(transpose(J)*J)
        ah[1,P] = â[1,1]; ah[2,P] = â[1,2]; ah[3,P] = â[2,2]
    end
    return ah
end

# ── BILINEAR straight-sided quad: the geometry model that matches real meshes ─
# A straight-sided general quad (trapezoid/kite/rhomboid) maps from the reference
# square by the BILINEAR (Q1) map  x(ξ,η) = Σ_i N_i(ξ,η) P_i,  with corners ordered
#   P1=(ξ,η)=(-1,-1),  P2=(+1,-1),  P3=(+1,+1),  P4=(-1,+1)
# and N_i = ¼(1±ξ)(1±η). This holds at ANY solution order: the high-order LGL
# nodes lie ON the straight edges, so the geometry is fixed by the 4 corners while
# â is still evaluated at all (k+1)² nodes (and the stiffness uses the order-k
# basis). J = ∂x/∂ξ = [∂x/∂ξ ∂x/∂η; ∂y/∂ξ ∂y/∂η] varies linearly inside the quad.
function bilinear_jacobian(P1, P2, P3, P4, ξ, η)
    dNξ1 = -(1-η)/4; dNξ2 =  (1-η)/4; dNξ3 =  (1+η)/4; dNξ4 = -(1+η)/4
    dNη1 = -(1-ξ)/4; dNη2 = -(1+ξ)/4; dNη3 =  (1+ξ)/4; dNη4 =  (1-ξ)/4
    xξ = dNξ1*P1[1] + dNξ2*P2[1] + dNξ3*P3[1] + dNξ4*P4[1]
    yξ = dNξ1*P1[2] + dNξ2*P2[2] + dNξ3*P3[2] + dNξ4*P4[2]
    xη = dNη1*P1[1] + dNη2*P2[1] + dNη3*P3[1] + dNη4*P4[1]
    yη = dNη1*P1[2] + dNη2*P2[2] + dNη3*P3[2] + dNη4*P4[2]
    return [xξ xη; yξ yη]
end

# Per-node â(ξ) field of a bilinear quad, evaluated at the LGL nodes ξnodes (the
# SAME nodes Jexpresso uses), in tensor flat order P=(i-1)*ngl+j (node (i,j) sits
# at (ξ_i, ξ_j)). â = a·det(J)·J⁻¹J⁻ᵀ — identical to el_ahat_nodes_from_metrics,
# so the sampled feature matches what inference builds from the real metrics.
# Returns (ah::Matrix(3,ngl²), ok::Bool) with ok=false if det(J)≤0 anywhere.
function ahat_field_from_quad(a, P1, P2, P3, P4, ξnodes, ngl::Int)
    ah = Matrix{Float64}(undef, 3, ngl*ngl)
    ok = true
    @inbounds for i = 1:ngl, j = 1:ngl
        P  = (i-1)*ngl + j
        J  = bilinear_jacobian(P1, P2, P3, P4, ξnodes[i], ξnodes[j])
        Je = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        Je <= 0 && (ok = false)
        Ji = inv(J)
        c  = a*Je
        ah[1,P] = c*(Ji[1,1]^2     + Ji[1,2]^2)        # â11
        ah[2,P] = c*(Ji[1,1]*Ji[2,1] + Ji[1,2]*Ji[2,2])# â12
        ah[3,P] = c*(Ji[2,1]^2     + Ji[2,2]^2)        # â22
    end
    return ah, ok
end

# Draw a random straight-sided quad: a base affine map (rotation × stretch × shear)
# applied to the reference corners, then an independent per-corner jitter to break
# parallelism (⇒ a genuine bilinear quad). Amplitude a ~ U(amin,amax). Returns
# (a, P1, P2, P3, P4).
function synthesize_random_quad(rng::AbstractRNG;
                                amin=1.0, amax=1.0,   # a=1: nuisance (Tie is amplitude-invariant), matches inference
                                smin=0.5, smax=2.0,
                                shear_max=0.5, jitter=0.25)
    a  = amin + (amax-amin)*rand(rng)
    θ  = 2π*rand(rng)
    s1 = smin + (smax-smin)*rand(rng)
    s2 = smin + (smax-smin)*rand(rng)
    γ  = shear_max*(2*rand(rng) - 1)
    B  = [cos(θ) -sin(θ); sin(θ) cos(θ)] * [s1 γ; 0.0 s2]
    refc = ((-1.0,-1.0), (1.0,-1.0), (1.0,1.0), (-1.0,1.0))
    P = map(refc) do r
        b = B*[r[1], r[2]]
        (b[1] + jitter*(2*rand(rng)-1), b[2] + jitter*(2*rand(rng)-1))
    end
    return a, P[1], P[2], P[3], P[4]
end

# ── T^{ie} = (A_{v^o,v^o})^{-1} A_{v^o,v^b} from a local stiffness matrix ─────
# Returns the (nint × nbdy) matrix in the convention used by the sampler
# (output_tensor stores vec of this; matches the sign used elsewhere in the EL
# code where the recorded output is +vec(A_vovo^{-1} A_vovb)).
function el_compute_Tie(A_loc::AbstractMatrix, vo::Vector{Int}, vb::Vector{Int})
    Avovo = A_loc[vo, vo]
    Avovb = A_loc[vo, vb]
    return Avovo \ Avovb
end

# ── Per-node â feature vector (length 3·ngl²) ────────────────────────────────
# Layout: for node P = 1..ngl², entries [3P-2, 3P-1, 3P] = (â11_P, â12_P, â22_P).
# `ahat_nodes` is (3, ngl²): row 1=â11, 2=â12, 3=â22 at each node.
function el_ahat_feature(ahat_nodes::AbstractMatrix)
    n = size(ahat_nodes, 2)
    f = Vector{Float64}(undef, 3n)
    @inbounds for P = 1:n
        f[3P-2] = ahat_nodes[1,P]
        f[3P-1] = ahat_nodes[2,P]
        f[3P]   = ahat_nodes[3,P]
    end
    return f
end

# Constant-â convenience: replicate (â11,â12,â22) across all ngl² nodes.
function el_ahat_feature_const(â11, â12, â22, ngl::Int)
    n = ngl*ngl
    ah = Matrix{Float64}(undef, 3, n)
    @inbounds for P = 1:n
        ah[1,P] = â11; ah[2,P] = â12; ah[3,P] = â22
    end
    return el_ahat_feature(ah)
end

# ── Per-node â from a REAL element's metrics (inference side) ─────────────────
# At node q the inverse Jacobian is J_K^{-1} = [dξdx dξdy; dηdx dηdy] and the
# Jacobian determinant is Je. With scalar physical diffusivity a_q the reference
# diffusivity is  â_q = Je_q · a_q · J_K^{-1} J_K^{-T}  (SPD). Returns (3, ngl²)
# in the P=(k-1)*ngl+l node order used throughout this module.
function el_ahat_nodes_from_metrics(dξdx, dξdy, dηdx, dηdy, Je, a_nodes, ngl::Int)
    n  = ngl*ngl
    ah = Matrix{Float64}(undef, 3, n)
    @inbounds for k = 1:ngl, l = 1:ngl
        P  = (k-1)*ngl + l
        ix = dξdx[k,l]; iy = dξdy[k,l]; nx = dηdx[k,l]; ny = dηdy[k,l]
        c  = Je[k,l]*a_nodes[k,l]
        ah[1,P] = c*(ix*ix + iy*iy)      # â11
        ah[2,P] = c*(ix*nx + iy*ny)      # â12
        ah[3,P] = c*(nx*nx + ny*ny)      # â22
    end
    return ah
end

# ── INFERENCE feature: per-element â from a real mesh's metrics ──────────────
# Scalar physical diffusivity a(x,y); default 1 (geometry-only â — the standard
# Laplacian amplitude). Redefine for a spatially-varying coefficient, e.g.
#   el_diffusivity(x,y) = 1.0 + (x > 5.0)   # cf. the PDF's a = 1 + 1(x>5).
el_diffusivity(x, y) = 1.0

# Fill `avisc` (nelem × 3·ngl²) with the per-node 2×2 SPD â feature for every
# element, in mesh.conn node order (matching el_nonconstant_sampling! and the
# elementLearning_infer! recovery). Built from the element metrics
# (dξdx,…,Je via J_K⁻¹) and the scalar diffusivity a evaluated at each node.
function el_avisc_nonconstant!(avisc::AbstractMatrix, mesh, metrics,
                               conn2ij, ngl::Int; afun = el_diffusivity)
    @inbounds for iel = 1:mesh.nelem
        for m = 1:ngl*ngl
            (i, j) = conn2ij[m]
            ix = metrics.dξdx[iel,i,j]; iy = metrics.dξdy[iel,i,j]
            nx = metrics.dηdx[iel,i,j]; ny = metrics.dηdy[iel,i,j]
            gid = mesh.connijk[iel, i, j, 1]
            a   = afun(mesh.x[gid], mesh.y[gid])
            c   = metrics.Je[iel,i,j] * a
            avisc[iel, 3m-2] = c*(ix*ix + iy*iy)   # â11
            avisc[iel, 3m-1] = c*(ix*nx + iy*ny)   # â12
            avisc[iel, 3m]   = c*(nx*nx + ny*ny)   # â22
        end
    end
    return avisc
end

# ── Sampling driver: generate (feature → Tie) training pairs (Option 1) ──────
# Produces, per sample, the per-node â feature (3·ngl²) and the flattened
# T^{ie}=(A_vovo)⁻¹A_vovb, BOTH in mesh.conn node ordering so the trained model
# is directly usable by elementLearning_infer! at inference time.
#
#   conn2ij      : conn-local index → tensor (i,j) map (from el_conn_to_ij)
#   elnbdypoints : number of boundary nodes per element (conn cols 1:elnbdypoints)
#   shape        : element-shape distribution used to draw â (see below)
#   ξnodes       : 1-D LGL nodes (length ngl); required for shape == :quad
#
# shape options:
#   :affine → random affine Jacobian (PARALLELOGRAM): â constant within the
#             element. Fast (uses the precomputed K-components).
#   :quad   → random straight-sided BILINEAR quad (4 corners): â varies within
#             the element exactly as a real mesh quad does. ★ recommended to
#             match straight-sided meshes (any solution order — see header).
#   :warp   → smooth synthetic per-node warp (legacy ξ-dependent approximation).
#
# Feature layout (per node, conn order m=1..ngl²): [â11_m, â12_m, â22_m].
# Output  layout: vec(T^{ie}) with T^{ie} of size (nvo × nvb) in conn order.
function el_nonconstant_sampling!(bufferin::Vector{Vector{Float64}},
                                  bufferout::Vector{Vector{Float64}},
                                  ψ, dψ, ω, ngl::Int, Nsamp::Int;
                                  conn2ij,
                                  elnbdypoints::Int,
                                  shape::Symbol = :affine,
                                  ξnodes::AbstractVector = Float64[],
                                  rng::AbstractRNG = Random.default_rng(),
                                  verbose::Bool = true)
    nelpoints = ngl*ngl
    perm = el_conn_perm(conn2ij, ngl)            # conn-local m → tensor flat P
    vb   = collect(1:elnbdypoints)               # conn order: boundary first
    vo   = collect(elnbdypoints+1:nelpoints)     #             interior last
    K11, K12, K22 = el_precompute_reference_components(ψ, dψ, ω, ngl)

    if shape == :quad && length(ξnodes) != ngl
        error("el_nonconstant_sampling!: shape=:quad needs ξnodes of length ngl=$ngl (got $(length(ξnodes)))")
    end

    for isamp = 1:Nsamp
        # 1) draw â in tensor order according to the chosen element-shape model
        â11 = â12 = â22 = 0.0
        ah_tensor = Matrix{Float64}(undef, 0, 0)
        if shape == :quad
            local ok
            for _try = 1:100                      # resample until det(J) > 0 everywhere
                a, P1, P2, P3, P4 = synthesize_random_quad(rng)
                ah_tensor, ok = ahat_field_from_quad(a, P1, P2, P3, P4, ξnodes, ngl)
                ok && break
            end
            A_loc = el_assemble_local_stiffness_field(ψ, dψ, ω, ah_tensor, ngl)
        elseif shape == :warp
            ah_tensor = synthesize_random_ahat_field(rng, ngl)
            A_loc     = el_assemble_local_stiffness_field(ψ, dψ, ω, ah_tensor, ngl)
        else  # :affine
            â11, â12, â22, _J = synthesize_random_ahat(rng)
            A_loc = el_assemble_local_stiffness(K11, K12, K22, â11, â12, â22)
        end

        # 2) reorder tensor → conn, partition, and condense  T^{ie}=Avovo⁻¹Avovb
        A_conn = A_loc[perm, perm]
        Tie    = A_conn[vo, vo] \ A_conn[vo, vb]

        # 3) build the per-node feature in conn order (node m ↔ tensor perm[m])
        feat = Vector{Float64}(undef, 3*nelpoints)
        @inbounds for m = 1:nelpoints
            if shape == :affine
                feat[3m-2] = â11; feat[3m-1] = â12; feat[3m] = â22
            else
                P = perm[m]
                feat[3m-2] = ah_tensor[1,P]; feat[3m-1] = ah_tensor[2,P]; feat[3m] = ah_tensor[3,P]
            end
        end

        push!(bufferin,  feat)
        push!(bufferout, vec(Tie))
        if verbose && (isamp % max(1, Nsamp ÷ 20) == 0)
            println(" # EL non-constant sampling ($(shape)): $isamp / $Nsamp")
        end
    end
    return length(vo), length(vb)
end
