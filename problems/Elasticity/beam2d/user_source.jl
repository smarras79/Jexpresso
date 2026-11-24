#---------------------------------------------------------------------------------
# S = (0, 0, 0, 0, 0, vx, vy)ᵀ
#
# The elastodynamic rows 1-5 have NO source — every term is a divergence. The
# load is not a body force here: it is a prescribed traction on the top
# surface, applied in user_bc.jl, which is also what lets it be ramped in time
# (user_source! is not handed t; user_bc_dirichlet! is). A gravity-like body
# force would go in S[2] as -ρg.
#
# S[6] and S[7] are the displacement recovery, u̇ₓ = vx and u̇_y = vy.
#
# This is the 2D source path (rhs.jl, NSD_2D), which weights S by ω·Jac exactly
# as it weights the flux — 2D and 3D always did. It is the 1D path that used to
# leave the Jacobian off; see the note above _expansion_inviscid!(..., NSD_1D,
# ContGal) in src/kernel/operators/rhs.jl.
#---------------------------------------------------------------------------------
function user_source!(S, q, qe, npoin::Int64, ::CL, ::TOTAL;
                      neqs=7, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    P = elastic_properties()

    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0
    S[5] = 0.0
    S[6] = q[1]/P.ρ     # u̇ₓ = vx
    S[7] = q[2]/P.ρ     # u̇_y = vy
end

#
# The PERT / NCL spellings exist only so that this file's `include` really does
# overwrite whatever a previously-run case left behind at those signatures (see
# the long note in problems/AdvDiff/kopriva/user_source.jl).
#
function user_source!(S, q, qe, npoin::Int64, ::CL, ::PERT;
                      neqs=7, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source!(S, q, qe, npoin::Int64, ::NCL, ::AbstractPert;
                      neqs=7, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    P = elastic_properties()
    return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(q[1]/P.ρ), T(q[2]/P.ρ)
end
