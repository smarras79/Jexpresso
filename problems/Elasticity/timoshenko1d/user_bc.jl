#---------------------------------------------------------------------------------
# Simply supported (pinned-pinned) beam.
#
# A pin holds the beam down but lets it rotate:
#
#     w = 0   ⟹  v = ẇ = 0        ⟹  U₁ = ρA v = 0
#     M = 0   ⟹  χ = M/EI = 0     ⟹  U₄ = χ     = 0
#
# ω and γ are left free — they are outgoing at both ends.
#
# That is exactly one condition per boundary per characteristic pair: the flux
# Jacobian splits into the (U₁,U₃) pair travelling at ±√(κₛG/ρ) and the (U₂,U₄)
# pair travelling at ±√(E/ρ), each of which has a single incoming characteristic
# at each end. Prescribing U₁ (first pair) and U₄ (second pair) is well posed.
#
# MECHANICS OF THE HOOK: the 1D driver prefills qbdy with a sentinel and only
# copies back the entries this function actually writes (see
# build_custom_bcs_dirichlet!(::NSD_1D, ...) in
# src/kernel/boundaryconditions/BCs.jl). Leaving qbdy[2] and qbdy[3] untouched
# is therefore how "no condition on ω and γ" is spelled — do not zero them.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)

    qbdy[1] = 0.0    # ρA v : the support holds w = 0
    qbdy[4] = 0.0    # χ    : a pin carries no bending moment

    return qbdy
end

function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::PERT)

    qbdy[1] = 0.0
    qbdy[4] = 0.0

    return qbdy
end

function user_bc_neumann(q::AbstractArray, gradq, coords, t, inputs)
    flux = zeros(size(q, 2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, lpert)
    T = eltype(q)
    return T(0.0), T(q[2]), T(q[3]), T(0.0)
end
