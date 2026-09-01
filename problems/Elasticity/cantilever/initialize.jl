function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the 1D Timoshenko beam (cantilever) ......... ")

    #-----------------------------------------------------------------------------
    # Conserved variables, in the order used by user_flux! / user_source!:
    #
    #     U = (ρA v, ρI ω, γ, χ, w, φ)ᵀ
    #
    # Rows 1-4 are the Timoshenko conservation laws; rows 5-6 are the passive
    # displacement/rotation recovery (zero flux, ẇ = v and φ̇ = ω in the source)
    # that turns the run into a picture of the beam actually bending.
    #
    # qoutvars are what HDF5/VTK write: the same state with the momentum
    # densities divided through by their inertia. The 1D PNG plotter plots the
    # CONSERVED set (qvars).
    #-----------------------------------------------------------------------------
    qvars    = ["ρAv", "ρIω", "γ", "χ", "w", "φ"]
    qoutvars = ["v",   "ω",   "γ", "χ", "w", "φ"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    #-----------------------------------------------------------------------------
    # Released from the static shape a uniform load q₀ holds the beam in, at
    # rest. The load is then removed (user_source! sets S[1] = 0), so the beam
    # swings about the undeformed configuration between ±w_static.
    #
    # Why release rather than suddenly APPLY the load: applying q₀ to a beam at
    # rest makes v grow uniformly everywhere except at the clamp, where the wall
    # holds v = 0. That is a step in x propagating in from the wall at the shear
    # speed, and a non-dissipative spectral element scheme rings on it. The
    # static shape, by contrast, is smooth, and it satisfies every boundary
    # condition of this case exactly — γ(L) = χ(L) = 0 at the free tip,
    # v = ω = w = φ = 0 at the clamp — so the initial data is compatible and the
    # solution stays smooth for the whole run.
    #-----------------------------------------------------------------------------
    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]
            x  = mesh.coords[1, ip]

            s = timoshenko_static_cantilever(x)

            q.qn[ip, 1] = 0.0        # ρA v : at rest
            q.qn[ip, 2] = 0.0        # ρI ω : at rest
            q.qn[ip, 3] = s.γ        # γ = Q/κₛGA
            q.qn[ip, 4] = s.χ        # χ = M/EI
            q.qn[ip, 5] = s.w        # w
            q.qn[ip, 6] = s.φ        # φ

            # No background state to subtract: the reference configuration of
            # this problem is the undeformed beam at rest.
            for ieq = 1:6
                q.qe[ip, ieq] = 0.0
            end
        end
    end

    p = timoshenko_properties()
    tip = timoshenko_static_cantilever(p.L)
    println("   beam: L = ", p.L, ", h = ", p.h, ", L/h = ", p.L/p.h)
    println("   wave speeds: shear √(κₛG/ρ) = ", p.cs, ", extensional √(E/ρ) = ", p.ce)
    println("   uniform load q₀ = ", p.q0,
            " → static tip deflection w(L) = ", tip.w,
            " (bending ", p.q0*p.L^4/(8.0*p.EI),
            " + shear ", p.q0*p.L^2/(2.0*p.κGA), ")")
    println(" Initialize fields for the 1D Timoshenko beam (cantilever) ......... DONE ")

    return q
end
