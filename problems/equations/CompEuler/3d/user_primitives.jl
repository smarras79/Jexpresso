function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
    uprimitive[5] = u[5]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = u[4]/(u[1]+qe[1])
    uprimitive[5] = (u[5]+qe[5])/(u[1]+qe[1])-qe[5]/qe[1]
end

function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T(u[4]/(u[1]+qe[1])), T((u[5]+qe[5])/(u[1]+qe[1]) - qe[5]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1]), T(u[5]/u[1])
    end
end


function user_uout!(uout, u, qe, ::TOTAL)

    #
    # IMPORTANT NOTICE:
    # uout has length equivalent to 'qoutvars' defined in your initialize.jl
    # If you have not defined qoutvars in initialize.jl, then uout will only
    # be allocated to the length of u by default.
    #
    uout[1] = u[1]
    uout[2] = u[2]/u[1]
    uout[3] = u[3]/u[1]
    uout[4] = u[4]/u[1]
    uout[5] = u[5]/u[1]
    
end


function user_uout!(uout, u, qe, ::PERT)

    uout[1] = u[1]-qe[1]
    uout[2] = u[2]/u[1]
    uout[3] = u[3]/u[1]
    uout[4] = u[4]/u[1]
    uout[5] = u[5]/u[1]
    
end
