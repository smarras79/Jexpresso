#---------------------------------------------------------------------------------
# Boundary conditions for the free-stream preservation test.
#
# TWO CONSTANTS ARE THE EXPERIMENT.  Flip one, rerun, compare.  They are
# deck constants and not environment switches so that a run is reproducible
# from the file alone.
#
#   FSP_OUTFLOW
#     :nothing     the PRODUCTION condition of rampCaoEtAl2021 -- supersonic
#                  outflow, impose nothing, let the interior convect out.
#                  This is the default and it is the hypothesis under test.
#     :freestream  the CONTROL -- prescribe the free stream on the outflow
#                  too.  Now EVERY boundary node is overwritten every stage,
#                  so a missing boundary term cannot show anywhere.
#
#   FSP_WALL
#     :freestream  the default.  The wall is REMOVED -- the bottom boundary
#                  is prescribed free stream like every other one, so there
#                  is no boundary layer, no no-slip discontinuity and no
#                  isothermal wall to confound the measurement.
#     :noslip      the production wall, for a second experiment once the
#                  first has been read.
#
# READING THE PAIR.  If :nothing shows a growing dp and :freestream shows
# machine zero, the defect is in what the scheme does at a boundary node
# where nothing is imposed -- i.e. a missing or wrong boundary term in the
# inviscid assembly -- and it is a BUG, not a stabilisation shortfall.  If
# BOTH show a growing dp, the defect is in the volume operator or the
# metrics and the outflow is innocent.  If NEITHER does, free-stream
# preservation is fine and the Mach-7 failures have another cause.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# THE TWO BOUNDARY CONSTANTS.  Exactly one of each pair is uncommented.
# The deck rungs live in user_inputs.jl; these are the two that cannot,
# because user_bc_dirichlet! is not handed `inputs`.
#---------------------------------------------------------------------------------

#
# OUTFLOW.  ALREADY ANSWERED -- leave it alone unless you are re-checking.
#
# RUNG 0 ran with :nothing, the production condition, and dp_rel came back
# at +/-5.8e-13 on the outflow plane like everywhere else.  The "impose
# nothing" outflow manufactures NOTHING, so the :freestream control below
# is no longer needed; it is kept only so the comparison can be redone.
#
const FSP_OUTFLOW = :nothing        # production, and the one that passed
#const FSP_OUTFLOW = :freestream    # control: prescribe the stream here too

#
# WALL.  THIS IS RUNG 2.  Comment the first line, uncomment the second.
#
# :freestream REMOVES the wall -- the plate and the ramp are prescribed
# free stream like every other boundary, so the uniform stream is an exact
# solution of the whole problem and dp_rel measures the discretisation and
# nothing else.  This is what RUNG 0 and RUNG 1 need.
#
# :noslip restores the production wall.  From that point on the uniform
# stream is NO LONGER an exact solution -- a 1725 m/s stream standing on a
# no-slip isothermal wall is a genuine discontinuity and dp_rel WILL become
# large.  That is not a failure.  What is under test from RUNG 2 on is not
# the magnitude of dp_rel but its LOCATION: a real disturbance is bounded
# by the fastest signal speed, u + c = 1949 m/s, which crosses the 60 mm to
# the top boundary in 3.1e-5 s, i.e. step 31,000.  Anything appearing at
# the top of the domain inside 500 steps did not travel there.
#
const FSP_WALL    = :freestream     # RUNG 0, RUNG 1  -- wall removed
#const FSP_WALL   = :noslip         # RUNG 2 and up   -- production wall

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "outflow" && FSP_OUTFLOW === :nothing
        # Impose nothing, exactly as the production deck does.
        return nothing
    end

    if tag == "wall" && FSP_WALL === :noslip
        PhysConst = PhysicalConst{Float64}()
        cv        = PhysConst.Rair/PhysConst.γm1
        ρ         = q[1]
        qbdy[2]   = 0.0
        qbdy[3]   = 0.0
        qbdy[4]   = ρ*cv*293.0
        return nothing
    end

    # Everything else -- inflow, top, and whichever of wall/outflow the
    # constants above did not claim -- is the free stream.
    ρ∞, u∞, v∞, p∞, T∞, ρE∞ = fsp_freestream()
    qbdy[1] = ρ∞
    qbdy[2] = ρ∞*u∞
    qbdy[3] = ρ∞*v∞
    qbdy[4] = ρE∞

    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,
                            qe, ::PERT)
    # This case runs TOTAL(); the method exists so the dispatch is complete.
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
    return T(qbdy[1]), T(qbdy[2]), T(qbdy[3]), T(qbdy[4])
end
