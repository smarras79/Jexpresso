#---------------------------------------------------------------------------------
# S = (0, (c³ - c - μ)/τ)ᵀ
#
#   S[1]  none: ∂ₜc = ∂ₓ²μ is a pure conservation law, so ∫c dx is conserved.
#   S[2]  the undifferentiated part of τ ∂ₜμ = c³ - c - ε²∂ₓ²c - μ. The
#         -ε²∂ₓ²c is the viscous term of equation 2 (user_inputs.jl).
#
# The PERT / NCL spellings exist only so that this file's `include` really does
# overwrite whatever a previously-run case left behind at those signatures (see
# the note in problems/AdvDiff/kopriva/user_source.jl). This case runs with
# CL() and TOTAL().
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=2, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    τ = ch_parameters().τ
    c = q[1]
    μ = q[2]

    S[1] = 0.0
    S[2] = (c*c*c - c - μ)/τ
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=2, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::NCL, ::AbstractPert;
                      neqs=2, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source_gpu(q, qe, x, PhysConst, xmax, xmin, lpert)
    T = eltype(q)
    τ = ch_parameters().τ

    return T(0.0), T((q[1]*q[1]*q[1] - q[1] - q[2])/τ)
end
