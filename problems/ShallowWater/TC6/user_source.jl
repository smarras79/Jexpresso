#---------------------------------------------------------------------------------
# TC6: f-plane Coriolis + manufactured-solution forcing (paper Eq. 42).
#   S_h  =         s_h(x,y,t)
#   S_hu = +f hv + s_hu(x,y,t)
#   S_hv = -f hu + s_hv(x,y,t)
# The forcing depends on time. user_source! has no time argument, so the current
# right-hand-side (RK stage) time is read from RHS_TIME[] (src/kernel/operators/rhs.jl).
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    s_h, s_hu, s_hv = tc6_sources(x, y, RHS_TIME[])
    S[1] = s_h
    S[2] =  TC6_F * q[3] + s_hu
    S[3] = -TC6_F * q[2] + s_hv
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    error(" # ERROR TC6: only SOL_VARS_TYPE => TOTAL() is implemented.")
end
