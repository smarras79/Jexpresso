function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the 2D simply supported beam ................ ")

    #-----------------------------------------------------------------------------
    # U = (ρvx, ρvy, σxx, σyy, σxy, ux, uy)ᵀ
    #
    # Rows 1-5 are the elastodynamic conservation laws; rows 6-7 are the
    # passive displacement recovery (zero flux, u̇ = v in the source) that turns
    # the output into a picture of the beam actually sagging.
    #
    # NAMING, because the output does not warn you: vx/vy are the particle
    # VELOCITY, ux/uy are the DISPLACEMENT. The sag you want to look at is uy.
    #-----------------------------------------------------------------------------
    # ASCII ONLY. These strings become VTK/HDF5 array names, and the Greek
    # σ does not survive the round trip — VisIt shows a field written as
    # "σyy" under the name "yy". Greek stays in the comments, where it reads
    # better and cannot break anything.
    qvars    = ["rho_vx", "rho_vy", "sxx", "syy", "sxy", "ux", "uy"]
    qoutvars = ["vx",     "vy",     "sxx", "syy", "sxy", "ux", "uy"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    #-----------------------------------------------------------------------------
    # EVERYTHING STARTS AT ZERO: undeformed, unstressed, at rest.
    #
    # That is what makes this case well behaved. Every boundary condition is
    # satisfied exactly at t = 0 — σ ≡ 0 leaves the top and bottom surfaces
    # traction-free, and the load is ramped up from zero rather than switched
    # on, so the loaded surface is consistent too. There is no incompatibility
    # anywhere for the scheme to ring on.
    #
    # The whole problem is then driven by the boundary: the load grows on the
    # top surface, the beam sags, overshoots the static deflection, and rings
    # about it (undamped apart from the small artificial viscosity in the deck).
    #-----------------------------------------------------------------------------
    for ip = 1:mesh.npoin
        for ieq = 1:7
            q.qn[ip, ieq] = 0.0
            q.qe[ip, ieq] = 0.0
        end
    end

    P = elastic_properties()
    w = beam_static_midspan()
    println("   material: ρ = ", P.ρ, ", E = ", P.E, ", ν = ", P.ν,
            "  →  λ = ", P.λ, ", μ = ", P.μ)
    println("   wave speeds: cp = ", P.cp, ", cs = ", P.cs)
    println("   beam: span L = ", P.L, ", depth h = ", P.h, ", L/h = ", P.L/P.h)
    println("   load: total ", P.Ptot, " over a patch of width ", P.aload,
            " at midspan, ramped over t = ", P.tramp)
    println("   expected static midspan sag (beam theory): ", w.total,
            "  = bending ", w.bending, " + shear ", w.shear,
            "  (shear is ", round(100*w.shear/w.total, digits=1), "% — this is why L/h = 5)")
    println("   the beam should OVERSHOOT to about twice that before ringing back")
    println(" Initialize fields for the 2D simply supported beam ................ DONE ")

    return q
end
