#---------------------------------------------------------------------------------
# Primitives and the output map for the free-stream preservation test.
#
# user_primitives! is the standard total-energy set; it is never exercised
# here because :lvisc => false, but the methods must exist.
#
# user_uout! IS the measurement.  It writes the DEPARTURE from the free
# stream, not the state, because the state is 760 Pa everywhere and a plot
# of it says nothing:
#
#   dp     = p - p_inf              [Pa]      the signed error
#   dp_rel = (p - p_inf)/p_inf      [-]       the number to quote
#
# In ParaView, colour by dp_rel and read its RANGE off the Information tab.
# That range, at the last output time, is the whole result of this case.
# 1e-14 is machine zero and means the free stream is preserved.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    ρ  = u[1]; ρu = u[2]; ρv = u[3]; ρE = u[4]
    uprimitive[1] = ρ
    uprimitive[2] = ρu/ρ
    uprimitive[3] = ρv/ρ
    uprimitive[4] = ρE/ρ - 0.5*(ρu*ρu + ρv*ρv)/(ρ*ρ)
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    ρ  = u[1] + qe[1]; ρu = u[2] + qe[2]; ρv = u[3] + qe[3]; ρE = u[4] + qe[4]
    uprimitive[1] = ρ
    uprimitive[2] = ρu/ρ
    uprimitive[3] = ρv/ρ
    uprimitive[4] = ρE/ρ - 0.5*(ρu*ρu + ρv*ρv)/(ρ*ρ)
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    ρ  = u[1]; ρu = u[2]; ρv = u[3]; ρE = u[4]
    return T(ρ), T(ρu/ρ), T(ρv/ρ), T(ρE/ρ - T(0.5)*(ρu*ρu + ρv*ρv)/(ρ*ρ))
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)

    PhysConst = PhysicalConst{Float64}()
    _, _, _, p∞, _, _ = fsp_freestream()

    ρ  = u[1]; ρu = u[2]; ρv = u[3]; ρE = u[4]

    p = PhysConst.γm1*(ρE - 0.5*(ρu*ρu + ρv*ρv)/ρ)

    uout[1] = ρ
    uout[2] = ρu/ρ
    uout[3] = ρv/ρ
    uout[4] = p - p∞            # dp     [Pa]
    uout[5] = (p - p∞)/p∞       # dp_rel [-]   <- THE measurement
end
