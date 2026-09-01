#---------------------------------------------------------------------------------
# Clamped-free (cantilever) plane-strain beam.
#
#   "clamped"  (x = 0)   the wall holds the material:
#                            u = v = 0,  and the recovered uₓ = u_y = 0.
#                        The stresses are FREE — the wall carries whatever
#                        traction it must, and that reaction is the answer.
#
#   "free"     (the other three sides)   nothing holds the surface, so the
#                        surface transmits nothing: the traction σ·n vanishes.
#                        On this axis-aligned mesh n is (±1,0) or (0,±1), so
#
#                            n = (±1, 0):  σxx = 0 and σxy = 0
#                            n = (0, ±1):  σyy = 0 and σxy = 0
#
#                        The velocities are free — the tip is where the beam
#                        moves most.
#
# COUNTING: the flux Jacobian in the direction n has eigenvalues ±cp, ±cs and
# 0, so each boundary has exactly two incoming characteristics and takes
# exactly two conditions. The clamp supplies u and v; the free surface supplies
# the two components of the traction. Rows 6-7 have zero characteristic speed
# and are ODEs at each node; pinning them at the wall is redundant with
# u = v = 0 there, and is done anyway so the recovered displacement cannot
# creep over tens of thousands of steps.
#
# The two corners at x = 0 are on a "clamped" edge AND a "free" edge, and the
# driver visits each edge separately, so those nodes get both sets. That is
# correct here: the two sets touch disjoint components of U.
#
# MECHANICS OF THE HOOK: the driver prefills qbdy with a sentinel and only
# copies back the entries this function actually writes (see
# build_custom_bcs_dirichlet!(::NSD_2D, ...) in
# src/kernel/boundaryconditions/BCs.jl). Leaving an entry untouched is how "no
# condition on this variable" is spelled — do not zero the free ones.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "clamped"
        qbdy[1] = 0.0        # ρu
        qbdy[2] = 0.0        # ρv
        qbdy[6] = 0.0        # uₓ
        qbdy[7] = 0.0        # u_y
    else  # "free": traction-free, σ·n = 0
        if abs(nx) > abs(ny)
            qbdy[3] = 0.0    # σxx
            qbdy[5] = 0.0    # σxy
        else
            qbdy[4] = 0.0    # σyy
            qbdy[5] = 0.0    # σxy
        end
    end

    return qbdy
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::PERT)
    return user_bc_dirichlet!(q, coords, t, tag, qbdy, nx, ny, qe, TOTAL())
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords,
                         t::AbstractFloat, tag::String, inputs)
    return zeros(size(q, 2), 1)
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords,
                         t::AbstractFloat, inputs)
    return zeros(size(q, 2), 1)
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    # Clamped end; the CPU dispatch above handles tag-based selection.
    return T(0.0), T(0.0), T(q[3]), T(q[4]), T(q[5]), T(0.0), T(0.0)
end
