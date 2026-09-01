#---------------------------------------------------------------------------------
# Clamped-free (cantilever) beam.
#
#   x = 0, CLAMPED   the wall holds displacement and rotation:
#                        w = 0, φ = 0   ⟹  v = 0, ω = 0
#                    so U₁ = ρA v = 0, U₂ = ρI ω = 0, and the recovered U₅ = w
#                    and U₆ = φ are pinned to 0 as well. γ and χ are free —
#                    the wall carries whatever shear force and moment it must.
#
#   x = L, FREE      nothing holds the tip, so the tip transmits nothing:
#                        Q = κₛGA γ = 0  ⟹  γ = 0
#                        M = EI χ    = 0  ⟹  χ = 0
#                    v, ω, w, φ are free — the tip is where the beam moves most.
#
# That is exactly one condition per boundary per characteristic pair: the flux
# Jacobian splits into the (U₁,U₃) pair travelling at ±√(κₛG/ρ) and the (U₂,U₄)
# pair travelling at ±√(E/ρ), each with a single incoming characteristic at each
# end. The clamp supplies U₁ and U₂ (one per pair), the free end supplies U₃ and
# U₄. Rows 5-6 have zero characteristic speed and are ODEs at each node; pinning
# them at the wall is redundant with v = ω = 0 there, and is done anyway so the
# recovered displacement cannot creep over 30 000 steps.
#
# MECHANICS OF THE HOOK: the 1D driver prefills qbdy with a sentinel and only
# copies back the entries this function actually writes (see
# build_custom_bcs_dirichlet!(::NSD_1D, ...) in
# src/kernel/boundaryconditions/BCs.jl). Leaving an entry untouched is how "no
# condition on this variable" is spelled — do not zero the free ones.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)

    if tag == "left"
        # clamped wall
        qbdy[1] = 0.0    # ρA v
        qbdy[2] = 0.0    # ρI ω
        qbdy[5] = 0.0    # w
        qbdy[6] = 0.0    # φ
    else  # "right"
        # free tip
        qbdy[3] = 0.0    # γ  ⟹  Q = 0
        qbdy[4] = 0.0    # χ  ⟹  M = 0
    end

    return qbdy
end

function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::PERT)
    return user_bc_dirichlet!(q, coords, t, tag, qbdy, qe, TOTAL())
end

function user_bc_neumann(q::AbstractArray, gradq, coords, t, inputs)
    flux = zeros(size(q, 2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, lpert)
    T = eltype(q)
    # Clamped end; the CPU dispatch above handles tag-based selection.
    return T(0.0), T(0.0), T(q[3]), T(q[4]), T(0.0), T(0.0)
end
