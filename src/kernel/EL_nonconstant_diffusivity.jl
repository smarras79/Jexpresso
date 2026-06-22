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
#  i.e. a 3·(k+1)^2 input vector. For an affine element â is ξ-independent, so
#  the per-node entries are identical; the per-node layout keeps the interface
#  ready for ξ-dependent (curved-element) â later.
#
#  Sampling source (Option 1): random element shapes are produced by
#  *synthesizing* random affine Jacobians J_K directly (no mesh file needed),
#  which induces the distribution of â matrices used to train the surrogate.
#
#  This whole path is gated by inputs[:lEL_nonconstant] (default false), so the
#  existing constant-amplitude pipeline and its trained model are untouched.
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

# ── Synthesize a random affine Jacobian → SPD reference diffusivity â ─────────
# Draw an element shape J_K and a scalar amplitude a ~ U(amin,amax); return the
# constant SPD tensor â = a·det(J)·(JᵀJ)^{-1} as its 3 unique entries, plus J.
# The shape model: rotation θ, anisotropic stretch (s1,s2), and a shear γ — a
# generic distribution of well-formed (det>0) elements.
function synthesize_random_ahat(rng::AbstractRNG;
                                amin=0.5, amax=1.5,
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

# ── Sampling driver: generate (feature → Tie) training pairs (Option 1) ──────
# Writes the per-node â feature (3·ngl²) and the flattened Tie (nint·nbdy) to
# the input/output buffers, one column per sample. Self-contained: needs only
# the reference basis (ψ, dψ, ω) and ngl.
function el_nonconstant_sampling!(bufferin::Vector{Vector{Float64}},
                                  bufferout::Vector{Vector{Float64}},
                                  ψ, dψ, ω, ngl::Int, Nsamp::Int;
                                  rng::AbstractRNG = Random.default_rng(),
                                  verbose::Bool = true)
    vo, vb = el_reference_node_partition(ngl)
    K11, K12, K22 = el_precompute_reference_components(ψ, dψ, ω, ngl)

    for isamp = 1:Nsamp
        â11, â12, â22, _J = synthesize_random_ahat(rng)
        A_loc = el_assemble_local_stiffness(K11, K12, K22, â11, â12, â22)
        Tie   = el_compute_Tie(A_loc, vo, vb)

        feat  = el_ahat_feature_const(â11, â12, â22, ngl)   # 3·ngl²
        out   = vec(Tie)                                    # nint·nbdy

        push!(bufferin,  feat)
        push!(bufferout, out)
        if verbose && (isamp % max(1, Nsamp ÷ 20) == 0)
            println(" # EL non-constant sampling: $isamp / $Nsamp")
        end
    end
    return length(vo), length(vb)
end
