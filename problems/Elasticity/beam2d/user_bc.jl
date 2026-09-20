#---------------------------------------------------------------------------------
# Simply supported beam with a load on the middle of its top surface.
#
# Three tagged boundaries, one condition per characteristic family at each. In
# the direction n the flux Jacobian splits into a P family — the normal
# velocity vn paired with the normal traction Tn — and an S family — the
# tangential velocity vt paired with the tangential traction Tt. Each carries
# one incoming characteristic, so each boundary takes exactly two conditions,
# and each condition must name one variable from ONE family:
#
#   "support"  the two end faces, n = (±1, 0).   HINGED:
#                  S family:  vy  = 0    the support holds the end at height 0
#                  P family:  σxx = 0    it does not push back axially
#              vx, σyy, σxy stay free — the end is free to slide axially and to
#              ROTATE, which is what makes it a hinge rather than a clamp. (A
#              clamp would instead pin both velocities and leave the stresses
#              free; that is the difference between this case and a cantilever.)
#
#   "top"      n = (0, 1).   LOADED SURFACE:
#                  P family:  σyy = -p(x,t)   the applied downward pressure
#                  S family:  σxy = 0         nothing drags it sideways
#
#   "bottom"   n = (0, -1).  FREE SURFACE:  σyy = 0, σxy = 0.
#
# Two rollers rather than a pin-and-roller leaves rigid-body horizontal
# translation unconstrained. That is deliberate: the loading is symmetric and
# purely vertical, so vx stays antisymmetric about midspan and the beam never
# drifts. Pinning vx at one end instead would put a spurious axial reaction
# into a problem that has no axial load.
#
# Rows 6-7 (ux, uy) have zero characteristic speed and are ODEs at each node;
# uy is pinned at the supports to match vy = 0 so the recovered deflection
# cannot creep there over tens of thousands of steps.
#
# MECHANICS OF THE HOOK: the driver prefills qbdy with a sentinel and only
# copies back the entries this function actually writes (see
# build_custom_bcs_dirichlet!(::NSD_2D, ...) in
# src/kernel/boundaryconditions/BCs.jl). Leaving an entry untouched is how "no
# condition on this variable" is spelled — do not zero the free ones.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "support"
        qbdy[2] = 0.0                                  # ρvy : held at height 0
        qbdy[3] = 0.0                                  # σxx : no axial reaction
        qbdy[7] = 0.0                                  # uy  : ditto, recovered
    elseif tag == "top"
        qbdy[4] = -beam_load_pressure(coords[1], t)    # σyy : applied pressure
        qbdy[5] = 0.0                                  # σxy : nothing drags it
    else  # "bottom": free surface
        qbdy[4] = 0.0                                  # σyy
        qbdy[5] = 0.0                                  # σxy
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
    # Support face; the CPU dispatch above handles tag-based selection.
    return T(q[1]), T(0.0), T(0.0), T(q[4]), T(q[5]), T(q[6]), T(0.0)
end
