#---------------------------------------------------------------------------------
# S = (q, κₛGA γ, -ω, 0, v, ω)ᵀ
#
#   S[1]  distributed transverse load q(x,t) per unit length. ZERO during this
#         run: the load only ever enters through the initial condition, which is
#         the static shape it holds the beam in (see initialize.jl). The run is
#         the release — the beam swings about the undeformed configuration.
#
#         A load that stays on would sit here. Note that user_source! is not
#         handed the time, so a time-dependent q(x,t) needs the driver to pass
#         `t` down; q(x) alone works as written.
#
#   S[2]  the shear force Q = κₛGA γ feeding the angular-momentum balance.
#   S[3]  the -ω of γ̇ = v_x - ω.
#   S[4]  none: χ̇ = ω_x is already a pure conservation law.
#   S[5]  ẇ = v      displacement recovery
#   S[6]  φ̇ = ω      rotation recovery
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=6, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=1.0)

    p = timoshenko_properties()

    v = q[1]/p.ρA
    ω = q[2]/p.ρI
    γ = q[3]

    S[1] = 0.0          # q(x,t): released beam, no load during the run
    S[2] = p.κGA*γ
    S[3] = -ω
    S[4] = 0.0
    S[5] = v
    S[6] = ω
end

#
# The PERT / NCL spellings exist only so that this file's `include` really does
# overwrite whatever a previously-run case left behind at those signatures (see
# the long note in problems/AdvDiff/kopriva/user_source.jl). This case is run
# with CL() and TOTAL().
#
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=6, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=1.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::NCL, ::AbstractPert;
                      neqs=6, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=1.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source_gpu(q, qe, x, PhysConst, xmax, xmin, lpert)
    T = eltype(q)
    p = timoshenko_properties()

    return T(0.0), T(p.κGA*q[3]), T(-q[2]/p.ρI), T(0.0), T(q[1]/p.ρA), T(q[2]/p.ρI)
end
