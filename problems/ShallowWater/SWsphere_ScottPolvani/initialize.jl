#---------------------------------------------------------------------------------
# SWsphere_ScottPolvani — initial condition.
#
# FORCED-DISSIPATIVE SHALLOW-WATER TURBULENCE ON THE SPHERE
#   Scott, R. K. & Polvani, L. M. (2007), J. Atmos. Sci. 64, 3158-3176.
#
# THE INITIAL STATE IS REST. Unlike the freely decaying problem of Cho &
# Polvani (1996), where an initial random flow is what evolves, the
# forced-dissipative system starts from nothing and is driven: the fluid is
# motionless with a uniform depth H, and every bit of motion that follows is
# put there by the random small-scale forcing (src/kernel/operators/
# sphere_forcing.jl) and organised by rotation and the deformation radius.
# The paper's Fig. 1 shows exactly this — the kinetic energy growing
# linearly from zero at the injection rate — and its equilibria depend on the
# forcing and dissipation, not on any initial condition.
#
# H is not a free number. The one parameter of the resting layer that matters
# is the POLAR DEFORMATION RADIUS
#
#   L_D = √(gH) / (2Ω) ,
#
# which the paper prescribes as a fraction of the planetary radius (0.025 for
# Jupiter and Saturn, 0.1 for Uranus/Neptune), so
#
#   H = (2Ω L_D a)² / g ,      φ_ref = gH = (2Ω L_D a)² .
#
# Note that g itself then drops out of the dynamics — φ = gh is what the
# equations carry — and only reappears when h is written out. Near the
# equator, where L_D diverges, the relevant scale is the EQUATORIAL
# deformation radius L_eq = √(a L_D); both are printed.
#
# STATE VECTOR — the conservative Cartesian form of Marras, Kopera & Giraldo
# (2015), Eq. (8), exactly as in SWsphere:
#
#   q = [φ, φu, φv, φw]ᵀ ,   φ = g h ,   (u,v,w) the FULL 3-D Cartesian velocity
#
# so at rest q = [φ_ref, 0, 0, 0] everywhere. q.qe holds the same rest state:
# it is the reference φ_ref that the radiative relaxation -ν_h(φ - φ_ref)
# pulls the height back towards (paper Eq. 15, applied to the perturbation
# height only), and the natural thing to difference against.
#
# ALSO DONE HERE: the planet's Ω and g are handed to user_source.jl and
# user_primitives.jl, which receive no `inputs`, through the typed module
# globals SP_OMEGA / SP_GRAVITY declared in user_source.jl.
#---------------------------------------------------------------------------------

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    if rank == 0
        println(" # ")
        println(" # INITIAL CONDITION: Scott & Polvani (2007) forced-dissipative turbulence — fluid at rest ...")
    end

    #---------------------------------------------------------------------------------
    # The planet, from the deck. sp_radius is what the deck asked for; the shell
    # actually built has mesh.radius (the reader rescales the .msh to
    # :sphere_radius), and the balance below uses THAT, so the two are checked
    # against each other rather than assumed equal.
    #---------------------------------------------------------------------------------
    planet = get(inputs, :sp_planet, :unknown)
    a      = Float64(mesh.radius)
    Ω      = Float64(get(inputs, :sp_Omega,   7.292e-5))
    g      = Float64(get(inputs, :sp_gravity, 9.80616))
    LD     = Float64(get(inputs, :sp_LD,      0.025))       # in planetary radii
    Ro_t   = Float64(get(inputs, :sp_Ro_target, NaN))
    T      = 2π/Ω

    a_deck = Float64(get(inputs, :sp_radius, a))
    if abs(a - a_deck)/a_deck > 1.0e-6
        error(string(" # ERROR initialize.jl: the shell radius is ", a, " m but the deck's planet has a = ",
                     a_deck, " m. Set :sphere_radius => a (user_inputs.jl does) so the grid is rescaled."))
    end

    # Ω and g for the pointwise callbacks (see user_source.jl)
    global SP_OMEGA   = Ω
    global SP_GRAVITY = g

    # H from the deformation radius; :sp_H is only a cross-check of the deck's
    # own arithmetic.
    H    = (2Ω*LD*a)^2/g
    φref = g*H

    #---------------------------------------------------------------------------------
    # Solution variables: q = [φ, φu, φv, φw]; q.qout carries the PRIMITIVE
    # fields for output, [h, u, v, w], with the velocity projected on the shell
    # by user_uout! at write time.
    #---------------------------------------------------------------------------------
    qvars    = ["phi", "phiu", "phiv", "phiw"]
    qoutvars = ["h", "u", "v", "w"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qoutvars)

    if (inputs[:backend] == CPU())
        for ip = 1:mesh.npoin
            q.qe[ip, 1] = φref
            q.qe[ip, 2] = 0.0
            q.qe[ip, 3] = 0.0
            q.qe[ip, 4] = 0.0

            q.qn[ip, 1] = φref
            q.qn[ip, 2] = 0.0
            q.qn[ip, 3] = 0.0
            q.qn[ip, 4] = 0.0

            q.qout[ip, 1] = H
            q.qout[ip, 2] = 0.0
            q.qout[ip, 3] = 0.0
            q.qout[ip, 4] = 0.0
        end
    else
        error(" # ERROR problems/ShallowWater/SWsphere_ScottPolvani/initialize.jl: only the CPU backend is implemented.")
    end

    #---------------------------------------------------------------------------------
    # Report the setup in the paper's terms, so a run can be placed on its
    # Tables 1-2 and Fig. 13 at a glance. The element count is reduced first: in
    # parallel mesh.nelem is this rank's share.
    #---------------------------------------------------------------------------------
    nelem_g = Int(mesh.nelem)
    if MPI.Comm_size(comm) > 1
        nelem_g = MPI.Allreduce(nelem_g, MPI.SUM, comm)
    end
    if rank == 0
        c    = sqrt(φref)                     # gravity-wave speed
        Leq  = sqrt(a*LD*a)                   # equatorial deformation radius, dimensional
        U_t  = 2a*Ω*Ro_t                      # target velocity scale
        nRh  = sqrt(1/(2Ro_t))                # Rhines wavenumber, paper Eq. (7) with α → 0
        α    = 2Ro_t/LD^2                     # paper's α = aU/(ΩL_D²) = 2Ro/L_D²
        nf   = Int(get(inputs, :forcing_nf, 0))
        ε0   = Float64(get(inputs, :forcing_epsilon, 0.0))

        @printf(" #   planet: %s   a = %.4e m   Ω = %.4e 1/s (T = %.3e s = %.3f h)   g = %.3f m/s²\n",
                String(planet), a, Ω, T, T/3600, g)
        println(" #   Scott & Polvani (2007): ", get(inputs, :sp_figure, "?"))
        @printf(" #   L_D/a = %.4f  →  H = %.2f m,  φ_ref = gH = %.4e m²/s²,  √(gH) = %.2f m/s\n",
                LD, H, φref, c)
        @printf(" #   L_D = %.3e m ,  L_eq = √(a L_D) = %.3e m (%.3f a)\n", LD*a, Leq, Leq/a)
        if isfinite(Ro_t)
            @printf(" #   target Ro = %.4f  →  U = 2aΩRo = %.2f m/s ;  n_Rh ~ √(1/2Ro) = %.1f ;  α = 2Ro/L_D² = %.2f %s\n",
                    Ro_t, U_t, nRh, α,
                    α >= 1 ? "(α ≥ 1: no Rossby-wave arrest at high latitude — jets confined to low latitudes, paper §5b)" :
                             "(α < 1: Rhines regime)")
        end
        @printf(" #   energy input ε₀ = %.4e m²/s³ = %.4e a²/T³, forced at n_f = %d\n",
                ε0, ε0*T^3/a^2, nf)
        # 6n² quads on the cubed sphere, n·nop LGL nodes per panel edge, four
        # panels around the equator
        npanel = round(Int, sqrt(nelem_g/6))
        naround = 4*npanel*Int(mesh.nop)
        @printf(" #   grid: %d elements per panel edge, nop = %d → %d nodes around the equator, %.1f per wavelength at n_f\n",
                npanel, Int(mesh.nop), naround, naround/max(nf, 1))
        if isfinite(Ro_t) && nf > 0 && nf < 1.5*nRh
            @warn string(" # initialize.jl: the forcing wavenumber n_f = ", nf,
                         " is not well above the Rhines wavenumber ", @sprintf("%.1f", nRh),
                         " of the target Ro: the jets would form at the forcing scale. ",
                         "Raise :forcing_nf (and refine the grid) or pick a planet with a larger Ro.")
        end
        println(" #   state: at rest, h = H everywhere; all motion comes from the forcing")
    end

    #---------------------------------------------------------------------------------
    # The state is trivially on the shell (no momentum at all), but the check is
    # kept so the case fails the same way SWsphere would if someone edits the
    # loop above into something that is not.
    #---------------------------------------------------------------------------------
    drift = sphere_normal_momentum(q.qn, mesh; ivar = 2)
    if MPI.Comm_size(comm) > 1
        drift = MPI.Allreduce(drift, MPI.MAX, comm)
    end
    drift == 0.0 ||
        error(string(" # ERROR initialize.jl: the initial momentum is not tangential to the shell: max|(φu)·x̂| = ", drift))

    if rank == 0
        println(" # INITIAL CONDITION: Scott & Polvani (2007) forced-dissipative turbulence ... DONE")
    end

    return q
end
