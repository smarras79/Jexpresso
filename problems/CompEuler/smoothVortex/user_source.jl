#---------------------------------------------------------------------------------
# No source terms: this is the homogeneous compressible Euler system, and
# the isentropic vortex is an exact solution of it with no forcing of any
# kind. Every slot is zero and user_inputs.jl sets :lsource => false.
#
# That absence is the point of this case. Its MHD counterpart
# (problems/MHD/smoothVortex) carries the GLM divergence-cleaning field ψ and
# its damping source, so its "exact" solution is only exact for the ideal
# system and the numerical ∇·B biases the measured error. Here there is no B,
# no ψ, no source — what is left is the discretization alone.
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=1.0, xmin=0.0, xmax=3.0)

    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0

end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=1.0, xmin=0.0, xmax=3.0)

    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0

end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    return T(0.0), T(0.0), T(0.0), T(0.0)
end
