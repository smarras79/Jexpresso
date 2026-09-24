#---------------------------------------------------------------------------------
# Boundary conditions for the Mach-7.7 compression ramp (Cao et al. 2021,
# Section 2.3).
#
# The `tag` argument is the gmsh physical-curve name carried by the edge,
# so these five names are exactly the groups declared in ramp15.geo:
#
#   "inflow"    left boundary, x = 0              free stream
#   "top"       upper boundary                    free stream
#   "outflow"   right boundary, x = x_end         extrapolation
#   "wall"      bottom (plate and ramp)           no slip, isothermal
#
# Free stream (inflow and top).  "The free stream condition is also
# prescribed at the upper computational boundary" -- Section 2.3.  At
# M = 7.7 every characteristic enters through the left boundary, so the
# whole conservative state is prescribed there; the top is far enough from
# the shock system (see ramp15.geo) that prescribing it there is a
# statement about the undisturbed stream, not about the solution.
#
# Outflow.  "An extrapolation condition is used at the outflow boundary".
# The flow leaves supersonically along the entire boundary except inside
# the boundary layer, so NOTHING is imposed: build_custom_bcs_dirichlet!
# pre-fills qbdy with the sentinel 4325789.0 and copies back only the slots
# this routine overwrites, so not touching qbdy on an "outflow" edge leaves
# the interior solution to convect out untouched.
#
# Wall.  "For the no-slip wall, isothermal conditions are specified with
# the wall temperature being 293 K".  Momentum is zeroed and the total
# energy is reset to its value at T_w for the density the interior
# solution brings to the wall,
#
#     rho E |_wall = rho cv T_w        (u = v = 0)
#
# Density itself is left alone: it is the one variable a wall does not
# constrain, and the continuity equation at the wall supplies it.
#
# NO UPSTREAM STRIP.  Section 2.3 puts 20 grid points in 1 mm of free
# stream ahead of the leading edge, and this deck did too until the strip
# turned out to be the single worst cell on the grid: conforming blocks
# force it to carry the wall-clustered dy = 8e-6 m, while its lower
# boundary is a symmetry line, so it has no boundary layer and runs the
# full 1726 m/s at the finest wall-normal spacing in the domain.  When
# DynSGS saturates there the viscous CFL it implies is ~3e-10 s.  Measured:
# with the strip, the case died at step 304 at exactly that node.  The
# strip is gone and the inflow now sits on the leading edge.
#
# THE LEADING EDGE (0,0) is the node two conditions claim: it opens the
# "wall" edge and closes the "inflow" edge.  The kernel walks the boundary
# edge by edge and writes each result straight into uaux, so the node is
# constrained twice and the later edge wins, whichever that is.  The wall
# condition is the physical one -- the leading edge is the first point of
# the plate, where the boundary layer starts and where the whole
# viscous-interaction pressure gradient of figure 2(b) originates -- so the
# inflow branch skips that node and lets the wall have it, deterministically.
#---------------------------------------------------------------------------------

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "inflow" || tag == "top"
        # The leading edge belongs to the wall, not to the inflow (see above).
        if tag == "inflow" && coords[2] < 1.0e-12
            return nothing
        end
        ρ∞, u∞, v∞, p∞, T∞, ρE∞ = ramp_freestream()
        qbdy[1] = ρ∞
        qbdy[2] = ρ∞*u∞
        qbdy[3] = ρ∞*v∞
        qbdy[4] = ρE∞

    elseif tag == "outflow"
        # Supersonic outflow: impose nothing.

    else
        # "wall": no slip, isothermal at T_w = 293 K.
        PhysConst = PhysicalConst{Float64}()
        cv        = PhysConst.Rair/PhysConst.γm1     # EOS-consistent cv

        ρ       = q[1]
        qbdy[2] = 0.0
        qbdy[3] = 0.0
        qbdy[4] = ρ*cv*ramp_Twall()
    end

    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,
                            qe, ::PERT)
    #
    # This case runs in TOTAL() mode (see user_inputs.jl).  The PERT()
    # method exists so the dispatch is complete; it applies free slip to
    # the perturbation momentum, which is what the other CompEuler cases
    # do.  It is NOT the isothermal wall -- do not switch :SOL_VARS_TYPE
    # on this case without writing that condition first.
    #
    qnl     = nx*(q[2] + qe[2]) + ny*(q[3] + qe[3])
    qbdy[2] = (q[2] + qe[2] - qnl*nx) - qe[2]
    qbdy[3] = (q[3] + qe[3] - qnl*ny) - qe[3]

end

#
# Neumann: nothing is imposed on any boundary.  Both overloads must exist
# (the solver calls both signatures).
#
# The viscous operator is assembled in weak form, so the boundary integral
# of the viscous flux is what these return.  Zero is right on every
# boundary here: the free-stream and outflow boundaries are far from any
# gradient worth carrying, and on the wall the momentum and energy slots
# are strong Dirichlet, which overrides whatever flux the weak form
# accumulated at those nodes.
#
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

#
# GPU path.  The KA kernels carry no boundary tag, so this can only apply
# ONE condition to every boundary node -- free slip, as in the other
# CompEuler cases.  That is not this case's wall, and it is not its
# free stream either, so run rampCaoEtAl2021 on the CPU backend.
#
function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    qnl = nx*q[2] + ny*q[3]
    u   = q[2] - qnl*nx
    v   = q[3] - qnl*ny
    return T(qbdy[1]), T(u), T(v), T(qbdy[4])
end
