#---------------------------------------------------------------------------------
# SWsphereDSGS — primitive variables and output variables.
#
# Solution (conservative) state, Marras/Kopera/Giraldo Eq. (8):
#
#   q    = [φ, φu, φv, φw]      φ = gh, (u,v,w) the CARTESIAN velocity
#
# Primitives — what the viscous operator δν∇²(φu) differentiates, and what the
# artificial-viscosity machinery expects to see. These stay CARTESIAN, because
# that is the frame the equations are written and differentiated in:
#
#   uprim = [φ, u, v, w]        Cartesian
#
# OUTPUT — what a human reads off the VTK file. These are NOT the Cartesian
# components: they are the velocity PROJECTED ONTO THE SHELL, in the local
# tangent basis,
#
#   uout  = [h, u_zonal, v_meridional, w_radial]
#
#     u_zonal = u·e_λ ,   e_λ = (-sin λ,        cos λ,       0    )   eastward
#     v_merid = u·e_φ ,   e_φ = (-sin φ cos λ, -sin φ sin λ, cos φ)   northward
#     w_rad   = u·n̂   ,   n̂  = x/r                                    outward
#
# Plotting the raw Cartesian uₓ over a sphere shows a dipole straddling the
# prime meridian — an artefact of the frame, not of the flow — where the zonal
# component shows the jet as the band it actually is. w_radial is the velocity
# form of the constraint and must be ~0 everywhere; it is written precisely so
# that "is the flow still on the shell?" is answerable from the output file.
#
# λ and φ arrive as keyword arguments from the caller (sphere_time_loop.jl and
# drivers.jl pass mesh.lon[ip] / mesh.lat[ip]); if they are absent the Cartesian
# components are written unchanged rather than silently producing nonsense.
#
# g is a literal here for the same reason ω is one in user_source.jl: these
# callbacks receive no `inputs`. Keep it in step with :galewsky_gravity if the
# deck ever changes it.
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

    g = 9.80616

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
