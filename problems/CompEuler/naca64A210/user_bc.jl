#---------------------------------------------------------------------------------
# Boundary conditions for the Mach-3 aerofoil.
#
# The `tag` argument is the gmsh physical-curve name carried by the edge, so
# these names are exactly the groups declared in naca64A210.geo:
#
#   "inflow"   left wall,  x = 0                   supersonic inflow
#   "outflow"  right wall, x = 3                   supersonic outflow
#   "bottom"   tunnel floor, y = -1                free slip
#   "top"      tunnel roof,  y = +1                free slip
#   "airfoil"  the section itself                  free slip
#
# Supersonic inflow: every characteristic enters the domain, so the whole
# conservative state is prescribed at the free stream.
#
# Supersonic outflow: every characteristic leaves the domain, so NOTHING is
# imposed. build_custom_bcs_dirichlet! pre-fills qbdy with the sentinel
# 4325789.0 and only copies back the slots this routine overwrites, so simply
# not touching qbdy on an "outflow" edge leaves the interior solution to convect
# out untouched.
#
# Solid walls: free slip (reflecting), i.e. the momentum vector is projected
# onto the wall tangent, u·n = 0, while ρ and ρE are left alone. This is an
# inviscid calculation, so it is the condition on the aerofoil as well as on the
# tunnel.
#
# THE TUNNEL WALLS ARE WALLS, not a far field. They are 1 m from an aerofoil of
# chord 0.6 m, so the bow shock reflects off them and comes back onto the wake —
# this is a wind-tunnel calculation, exactly as shock_circle is, and those
# reflections are part of the answer rather than an artefact to be explained
# away. It is also why the free stream is horizontal and the incidence is zero:
# a free stream at an angle would fight a free-slip wall that is forcing v = 0.
# Putting an aerofoil in free air needs a characteristic far-field condition,
# which is a different boundary condition, not a different mesh.
#
# WHY THE AEROFOIL WALL IS WHERE THE CURVED GEOMETRY PAYS. The slip projection
# uses (nx, ny), the OUTWARD NORMAL OF THE DISCRETE WALL, which
# build_metric_terms! takes from the tangent of the degree-:nop interpolant
# through the boundary nodes. On a straight-sided grid that interpolant is the
# chord, so the normal is the polygon's — wrong by O(h) at every node, constant
# along each edge, and discontinuous at every vertex. Reflecting momentum in a
# normal that jumps around the nose of an aerofoil in a Mach-3 stream generates
# vorticity at each of those jumps. With :exact_geometry the interpolant is the
# true section, and the normal is the aerofoil's.
#
# CORNER NODES. build_custom_bcs_dirichlet! walks the boundary edge by edge and
# writes each result straight back into uaux, so a node shared by two edges is
# constrained TWICE, the second projection acting on the state the first one
# already modified. Here that is harmless everywhere it happens:
#
#   (0,±1), tunnel corners at inflow   the free stream has v = 0, which is
#                                      exactly what the slip projection with
#                                      n = (0,∓1) leaves behind — they agree.
#   (3,±1), tunnel corners at outflow  outflow imposes nothing, so only the
#                                      wall acts.
#
# The aerofoil's own trailing edge is a convex corner and would be the exception
# — two projections there would zero both momentum components and plant a
# stagnation point in the middle of the recompression — but it is not a corner
# of the MESH in that sense: it is a single node shared by two boundary edges
# whose normals differ by the (small) trailing-edge angle, not by 90° as on the
# forward-facing step of ffs_step. The two projections nearly agree, and the
# residual is the physical stagnation the trailing edge really does produce.
#---------------------------------------------------------------------------------

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "inflow"
        ρ∞, u∞, v∞, p∞, ρE∞ = naca_freestream()
        qbdy[1] = ρ∞
        qbdy[2] = ρ∞*u∞
        qbdy[3] = ρ∞*v∞
        qbdy[4] = ρE∞

    elseif tag == "outflow"
        # Supersonic outflow: impose nothing.

    else
        # "bottom", "top", "airfoil": free slip. Remove the wall-normal
        # momentum component.
        qnl     = nx*q[2] + ny*q[3]
        qbdy[2] = q[2] - qnl*nx
        qbdy[3] = q[3] - qnl*ny
    end

    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,
                            qe, ::PERT)
    #
    # This case runs in TOTAL() mode (see user_inputs.jl). The PERT() method
    # exists so the dispatch is complete; it applies free slip to the
    # perturbation momentum, which is what the other CompEuler cases do.
    #
    qnl     = nx*(q[2] + qe[2]) + ny*(q[3] + qe[3])
    qbdy[2] = (q[2] + qe[2] - qnl*nx) - qe[2]
    qbdy[3] = (q[3] + qe[3] - qnl*ny) - qe[3]

end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords,
                         t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords,
                         t::AbstractFloat, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    qnl = nx*q[2] + ny*q[3]
    u   = q[2] - qnl*nx
    v   = q[3] - qnl*ny
    return T(qbdy[1]), T(u), T(v), T(qbdy[4])
end
