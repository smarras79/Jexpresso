#---------------------------------------------------------------------------------
# Boundary conditions for the Mach-3 forward-facing step.
#
# The `tag` argument is the gmsh physical-curve name carried by the edge, so
# these three names are exactly the groups declared in
# ffs_step_transfinite.geo:
#
#   "inflow"   left wall, x = 0            supersonic inflow
#   "outflow"  right wall, x = 3           supersonic outflow
#   "wall"     floor, step face, step top, tunnel roof     free slip
#
# Supersonic inflow: every characteristic enters the domain, so the whole
# conservative state is prescribed at the free stream.
#
# Supersonic outflow: every characteristic leaves the domain, so NOTHING is
# imposed. build_custom_bcs_dirichlet! pre-fills qbdy with the sentinel
# 4325789.0 and only copies back the slots this routine overwrites, so
# simply not touching qbdy on an "outflow" edge leaves the interior
# solution to convect out untouched.
#
# Solid walls: free slip (reflecting), i.e. the momentum vector is projected
# onto the wall tangent, u·n = 0, while ρ and ρE are left alone. This is the
# condition used for every wall of the tunnel — including the two faces of
# the step — in Woodward & Colella and in Nazarov & Hoffman Section 5.1.
#
# Corner nodes belong to two edges and are visited once per edge, last write
# winning. The two inflow corners (0,0) and (0,1) are the only nodes shared
# between an inflow edge and a wall edge, and there the two conditions agree:
# the free stream has v = 0, which is exactly what the slip projection with
# n = (0,∓1) leaves behind.
#---------------------------------------------------------------------------------

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "inflow"
        ρ∞, u∞, v∞, p∞, ρE∞ = ffs_freestream()
        qbdy[1] = ρ∞
        qbdy[2] = ρ∞*u∞
        qbdy[3] = ρ∞*v∞
        qbdy[4] = ρE∞

    elseif tag == "outflow"
        # Supersonic outflow: impose nothing.

    else
        # "wall": free slip. Remove the wall-normal momentum component.
        qnl     = nx*q[2] + ny*q[3]
        qbdy[2] = q[2] - qnl*nx
        qbdy[3] = q[3] - qnl*ny
    end

end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,
                            qe, ::PERT)
    #
    # This case runs in TOTAL() mode (see user_inputs.jl). The PERT()
    # method exists so the dispatch is complete; it applies free slip to
    # the perturbation momentum, which is what the other CompEuler cases do.
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
