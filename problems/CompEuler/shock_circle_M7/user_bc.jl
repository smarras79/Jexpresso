#---------------------------------------------------------------------------------
# Boundary conditions for the Mach-7 laminar cylinder.
#
# The `tag` argument is the gmsh physical-curve name, so these are exactly the
# groups of cylinder_M7.geo:
#
#   "inflow"    left edge,  x = 0        free stream (supersonic: all in)
#   "top"       y = +1                   free stream
#   "bottom"    y = -1                   free stream
#   "outflow"   right edge, x = 3        supersonic outflow: impose NOTHING
#   "cylinder"  the circle               NO SLIP, ISOTHERMAL at T_w = 300 K
#
# WHY top AND bottom CARRY THE FREE STREAM rather than a slip wall. The bow
# shock off a cylinder asymptotes to the Mach angle, asin(1/7) = 8.2°, so from
# the body it reaches |y| ≈ 0.2 + 2·tan(8.2°) ≈ 0.49 at the outflow — it never
# touches y = ±1. Prescribing the undisturbed free stream there is therefore
# exact, and unlike a slip wall it cannot reflect anything back into the shock
# layer.
#
# THE CYLINDER. No slip with an isothermal wall, the condition
# rampCaoEtAl2021 uses and the one that makes a surface heat flux meaningful:
#
#     u = v = 0,      ρE|wall = ρ cv T_w
#
# Density itself is left alone — it is the one variable a wall does not
# constrain, and continuity at the wall supplies it. cv is taken from the
# code's own PhysConst as Rair/(γ-1) so it is consistent with the equation of
# state the fluxes use; writing 718 here instead would put the wall energy and
# the pressure on slightly different gases.
#
# NOTE ON THE WALL NORMAL. There is no normal in the condition above — no slip
# needs none. What DOES need the geometry to be right is the viscous stress and
# the heat flux the solver computes at these nodes, and those are only as good
# as the wall the nodes lie on. That is why the deck carries
# :exact_geometry => Dict("cylinder" => (:circle, 1.0, 0.0, 0.2)): without it
# the high-order nodes sit on the CHORDS of a 64-sided polygon and the
# wall-normal temperature gradient — i.e. the heat flux — is wrong by O(h) at
# every node. On a curved wall that snap is not an optimisation, it is the
# default expectation.
#
# CORNERS. The four box corners are shared by two free-stream edges, which
# prescribe the same state, so the double application is harmless. The cylinder
# touches no other group.
#---------------------------------------------------------------------------------

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "inflow" || tag == "top" || tag == "bottom"
        ρ∞, u∞, v∞, p∞, T∞, ρE∞ = cyl_freestream()
        qbdy[1] = ρ∞
        qbdy[2] = ρ∞*u∞
        qbdy[3] = ρ∞*v∞
        qbdy[4] = ρE∞

    elseif tag == "outflow"
        # Supersonic outflow: impose nothing. build_custom_bcs_dirichlet!
        # pre-fills qbdy with a sentinel and copies back only what this
        # routine overwrites, so not touching it lets the solution convect out.

    else
        # "cylinder": no slip, isothermal at T_w.
        PhysConst = PhysicalConst{Float64}()
        cv        = PhysConst.Rair/PhysConst.γm1     # EOS-consistent cv

        ρ       = q[1]
        qbdy[2] = 0.0
        qbdy[3] = 0.0
        qbdy[4] = ρ*cv*cyl_Twall()
    end

    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,
                            qe, ::PERT)
    #
    # This case runs in TOTAL() mode (see user_inputs.jl). The PERT() method
    # exists so the dispatch is complete; it applies free slip to the
    # perturbation momentum, as the other CompEuler cases do. It is NOT the
    # isothermal wall — do not switch :SOL_VARS_TYPE on this case without
    # writing that condition first.
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
    return T(qbdy[1]), T(0.0), T(0.0), T(qbdy[4])
end
