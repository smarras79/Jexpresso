#---------------------------------------------------------------------------------
# S = (0, -w/τ, v/τ)ᵀ
#
#   S[1]  none: the u row is a pure conservation law.
#   S[2]  the -w of τ ∂ₜv = ∂ₓv - w.
#   S[3]  the +v of τ ∂ₜw = v - ∂ₓu.
#
# The (v, w) block of the source is a rotation with rate 1/τ (eigenvalues ±i/τ),
# so the explicit integrator needs a stability region that covers part of the
# imaginary axis; the low-storage RK4 in user_inputs.jl does.
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
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    τ = kdv_parameters().τ

    S[1] = 0.0
    S[2] = -q[3]/τ
    S[3] =  q[2]/τ
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::NCL, ::AbstractPert;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source_gpu(q, qe, x, PhysConst, xmax, xmin, lpert)
    T = eltype(q)
    τ = kdv_parameters().τ

    return T(0.0), T(-q[3]/τ), T(q[2]/τ)
end
