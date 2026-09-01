function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for 2D elastodynamics (periodic plane P-wave) ... ")

    #-----------------------------------------------------------------------------
    # U = (ρu, ρv, σxx, σyy, σxy)ᵀ — the pure velocity-stress system, no source
    # and no reconstructed extras. qoutvars divide the momenta through by ρ so
    # the VTK files carry the particle velocity itself.
    #-----------------------------------------------------------------------------
    # `vx`/`vy` are the particle VELOCITY; qoutvars divide the momenta through
    # by ρ so the VTK files carry it directly.
    # ASCII ONLY: these become VTK/HDF5 array names and the Greek σ does not
    # survive the round trip (VisIt shows "σyy" as "yy").
    qvars    = ["rho_vx", "rho_vy", "sxx", "syy", "sxy"]
    qoutvars = ["vx",     "vy",     "sxx", "syy", "sxy"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    for ip = 1:mesh.npoin
        x, y = mesh.x[ip], mesh.y[ip]
        s = elastic_plane_wave_state(x, y, 0.0)

        q.qn[ip, 1] = s.ρu
        q.qn[ip, 2] = s.ρv
        q.qn[ip, 3] = s.σxx
        q.qn[ip, 4] = s.σyy
        q.qn[ip, 5] = s.σxy

        # No background state: the reference configuration is the undeformed
        # body at rest.
        for ieq = 1:5
            q.qe[ip, ieq] = 0.0
        end
    end

    P = elastic_properties()
    W = elastic_plane_wave()
    println("   material: ρ = ", P.ρ, ", E = ", P.E, ", ν = ", P.ν,
            "  →  λ = ", P.λ, ", μ = ", P.μ)
    println("   wave speeds: cp = ", P.cp, ", cs = ", P.cs, "  (cp/cs = ", P.cp/P.cs, ")")
    println("   P-wave along n = (", W.nx, ", ", W.ny, "), |k| = ", W.k,
            ", ω = ", W.ω, ", period = ", 2π/W.ω)
    println(" Initialize fields for 2D elastodynamics (periodic plane P-wave) ... DONE ")

    return q
end
