#---------------------------------------------------------------------------------
# SWsphere_ScottPolvani_DynSGS — primitive variables and output variables.
#
# The SAME as problems/ShallowWater/SWsphere/user_primitives.jl, except that g
# is the planet's (the module global SP_GRAVITY declared in user_source.jl and
# set by initialize() from the deck) rather than the Earth's literal.
#
# Solution (conservative) state:   q     = [φ, φu, φv, φw]   φ = gh, Cartesian u
# Primitives (what ν∇ₛ² differentiates): uprim = [φ, u, v, w]  Cartesian
# OUTPUT (what a human reads off the VTK file): the velocity PROJECTED ONTO THE
# SHELL, in the local tangent basis,
#
#   uout  = [h, u_zonal, v_meridional, w_radial]
#
#     u_zonal = u·e_λ ,   e_λ = (-sin λ,        cos λ,       0    )   eastward
#     v_merid = u·e_φ ,   e_φ = (-sin φ cos λ, -sin φ sin λ, cos φ)   northward
#     w_rad   = u·n̂   ,   n̂  = x/r                                    outward
#
# w_radial is the velocity form of the tangency constraint and must be ~0.
# u_zonal is the field the paper is read in (its zonal mean against latitude
# is also written as a text file by the time loop, see :lzonal_mean).
#---------------------------------------------------------------------------------

function user_primitives!(u, qe, uprimitive, ::TOTAL)
    φ = u[1]
    uprimitive[1] = φ
    uprimitive[2] = u[2]/φ
    uprimitive[3] = u[3]/φ
    uprimitive[4] = u[4]/φ
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    user_primitives!(u, qe, uprimitive, TOTAL())
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    φ = u[1]
    return T(φ), T(u[2]/φ), T(u[3]/φ), T(u[4]/φ)
end


function user_uout!(ip, ::TOTAL, uout, u, qe; lon = nothing, lat = nothing, kwargs...)

    g = SP_GRAVITY

    φ  = u[1]
    ux = u[2]/φ
    uy = u[3]/φ
    uz = u[4]/φ

    uout[1] = φ/g                      # h [m]

    if lon === nothing || lat === nothing
        # no coordinates supplied: fall back to the Cartesian components rather
        # than pretending they are zonal/meridional
        uout[2] = ux
        uout[3] = uy
        uout[4] = uz
        return
    end

    sλ, cλ = sincos(lon)
    sφ, cφ = sincos(lat)

    uout[2] =  -ux*sλ    +  uy*cλ                    # u·e_λ  zonal      (eastward)
    uout[3] =  -ux*sφ*cλ -  uy*sφ*sλ + uz*cφ         # u·e_φ  meridional (northward)
    uout[4] =   ux*cφ*cλ +  uy*cφ*sλ + uz*sφ         # u·n̂    radial     (≈ 0)
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    user_uout!(ip, TOTAL(), uout, u, qe; kwargs...)
end
