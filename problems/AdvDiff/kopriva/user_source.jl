#
# NO SOURCE: pure advection-diffusion, one equation.
#
# WHY THE SAME BODY IS SPELLED OUT THREE TIMES rather than once against
# ::AbstractPert, which is what AdvDiff/2d does.
#
# src/run.jl loads a case by `include`ing its user_*.jl, and an include only
# REPLACES a method whose signature matches EXACTLY. It cannot remove a method
# the previously-run case defined at a different signature. So after, say,
#
#     Jexpresso.run_case("CompEuler", "theta")     # user_source!(…, ::CL, ::PERT)
#     Jexpresso.run_case("AdvDiff",   "kopriva")   # user_source!(…, ::CL, ::AbstractPert)
#
# BOTH methods are live, and since PERT <: AbstractPert, theta's is the more
# specific one and wins dispatch. kopriva runs with :SOL_VARS_TYPE => PERT() and
# neqs = 1, so theta's four-equation body ran off the end of the source vector:
#
#     BoundsError: attempt to access 1-element view(...) at index [2]
#       user_source! @ problems/CompEuler/theta/user_source.jl:37
#
# Covering ::TOTAL, ::PERT and (::NCL, ::AbstractPert) explicitly means every
# signature a sibling case might have used is one this file overwrites, so the
# include really does take effect and no leftover method can outrank these.
#
# (An earlier fix evicted the outgoing case's methods with Base.delete_method in
# run.jl instead. That is the tidier idea and it must not come back: deleting a
# hook that has been inlined into the precompiled RHS takes the session down
# with a StackOverflowError in type inference. See the note in src/run.jl.)
#
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0

end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0

end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::NCL, ::AbstractPert;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0

end

function user_source_gpu(q,qe,x,y,PhysConst, xmax, xmin, ymax, ymin,lpert)

    T = eltype(q)
    #
    # S(q(x)) = -ρg
    #
    return T(0.0)
end
