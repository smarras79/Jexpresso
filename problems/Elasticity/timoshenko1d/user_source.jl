#---------------------------------------------------------------------------------
# S = (q, κₛGA γ, -ω, 0)ᵀ
#
#   S[1]  distributed transverse load q(x,t) per unit length. Zero here: this
#         case is a free standing vibration released from rest.
#   S[2]  the shear force Q = κₛGA γ feeding the angular-momentum balance.
#   S[3]  the -ω of γ̇ = v_x - ω.
#   S[4]  none: χ̇ = ω_x is already a pure conservation law.
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=1.0)

    p = timoshenko_properties()

    ω = q[2]/p.ρI
    γ = q[3]

    S[1] = 0.0
    S[2] = p.κGA*γ
    S[3] = -ω
    S[4] = 0.0
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
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=1.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::NCL, ::AbstractPert;
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=1.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source_gpu(q, qe, x, PhysConst, xmax, xmin, lpert)
    T = eltype(q)
    p = timoshenko_properties()

    return T(0.0), T(p.κGA*q[3]), T(-q[2]/p.ρI), T(0.0)
end
