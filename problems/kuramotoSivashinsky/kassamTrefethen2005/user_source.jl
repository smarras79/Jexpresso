#---------------------------------------------------------------------------------
# S = (-w, -w/τ)ᵀ
#
#   S[1]  the -w of ∂ₜu + ∂ₓ(u²/2) = -w - ∂ₓ²w, i.e. -∂ₓ²u as τ → 0. The
#         -∂ₓ²w (→ -∂ₓ⁴u) is the viscous term of equation 1.
#   S[2]  the -w of τ ∂ₜw = ∂ₓ²u - w. The ∂ₓ²u is the viscous term of
#         equation 2.
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

    τ = ks_parameters().τ
    w = q[2]

    S[1] = -w
    S[2] = -w/τ
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
    τ = ks_parameters().τ

    return T(-q[2]), T(-q[2]/τ)
end
