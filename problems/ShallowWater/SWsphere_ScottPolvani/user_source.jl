#---------------------------------------------------------------------------------
# SWsphere_ScottPolvani — source terms of the shallow water equations on a
# spherical shell: Coriolis and the Lagrange multiplier.
#
# The SAME two terms as problems/ShallowWater/SWsphere/user_source.jl — see
# that file for the derivation of both — with one difference: the rotation
# rate is the PLANET'S, not the Earth's.
#
#   S_mom = -f (x × φu) + μ x ,     f = 2Ω z/r² ,     μ = -φ|u|²/r² .
#
# The forcing φu_F and the large-scale dissipation of Scott & Polvani (2007)
# are NOT here: user_source! is pointwise and has no access to a random field
# drawn once per step, so they live in src/kernel/operators/sphere_forcing.jl
# and enter the RHS after assembly (switched by :lsphere_forcing in the deck).
#
# WHERE Ω COMES FROM. user_source! receives no `inputs`, which is why SWsphere
# writes ω as a literal. Here the planet is a deck choice, so Ω is a TYPED
# MODULE GLOBAL, declared below with a placeholder and assigned by
# initialize() from :sp_Omega before the time loop starts. A typed global is
# read type-stably — the annotation is what keeps this hot loop free of
# dynamic dispatch; do not drop it — and unlike a `const` it can be reassigned
# when the deck switches planet in the same Julia session.
#---------------------------------------------------------------------------------

# Placeholders: Jupiter. initialize() overwrites both from the deck.
SP_OMEGA::Float64   = 1.7585e-4      # rotation rate [1/s]
SP_GRAVITY::Float64 = 24.79          # gravity [m/s²], used by user_primitives.jl


function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=4, x=0.0, y=0.0, z=0.0,
                      xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0, zmin=0.0, zmax=0.0)

    ω = SP_OMEGA

    φ  = q[1]
    φu = q[2]
    φv = q[3]
    φw = q[4]

    r2 = x*x + y*y + z*z

    #
    # Coriolis: -f (x × φu),  f = 2ω z/r²
    #
    f = 2.0*ω*z/r2

    cx = y*φw - z*φv
    cy = z*φu - x*φw
    cz = x*φv - y*φu

    #
    # Lagrange multiplier: μ = -φ|u|²/r², i.e. a centripetal force of magnitude
    # φ|u|²/r directed at the centre of the sphere.
    #
    u2 = (φu*φu + φv*φv + φw*φw)/(φ*φ)
    μ  = -φ*u2/r2

    S[1] = 0.0
    S[2] = -f*cx + μ*x
    S[3] = -f*cy + μ*y
    S[4] = -f*cz + μ*z
end


function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=4, x=0.0, y=0.0, z=0.0,
                      xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0, zmin=0.0, zmax=0.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, z=z,
                 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
end


function user_source_gpu(q, qe, x, y, z, PhysConst, lpert)
    T = eltype(q)

    ω  = T(SP_OMEGA)

    φ  = q[1]
    φu = q[2]
    φv = q[3]
    φw = q[4]

    r2 = x*x + y*y + z*z
    f  = T(2.0)*ω*z/r2

    cx = y*φw - z*φv
    cy = z*φu - x*φw
    cz = x*φv - y*φu

    u2 = (φu*φu + φv*φv + φw*φw)/(φ*φ)
    μ  = -φ*u2/r2

    return T(0.0), T(-f*cx + μ*x), T(-f*cy + μ*y), T(-f*cz + μ*z)
end
