function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the 1D Timoshenko beam (simply supported) ... ")

    #-----------------------------------------------------------------------------
    # Conserved variables, in the order used by user_flux! / user_source!:
    #
    #     U = (ρA v, ρI ω, γ, χ)ᵀ
    #
    # qoutvars are what HDF5/VTK write: the same state divided through by its
    # inertia, i.e. the physical velocity and angular rate. The 1D PNG plotter
    # plots the CONSERVED set (qvars) — that is the set user_analytic.jl is
    # asked for.
    #-----------------------------------------------------------------------------
    qvars    = ["ρAv", "ρIω", "γ", "χ"]
    qoutvars = ["v",   "ω",   "γ", "χ"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    m = timoshenko_mode()

    #-----------------------------------------------------------------------------
    # Released from rest at maximum deflection:
    #
    #     w = W sin(kx) cos(Ωt)        φ = Φ cos(kx) cos(Ωt)
    #
    # so at t = 0
    #
    #     v = ẇ = 0                    ω = φ̇ = 0
    #     γ = w_x - φ = (Wk - Φ) cos(kx)
    #     χ = φ_x     = -Φ k sin(kx)
    #
    # Only the flexural branch of the Timoshenko frequency equation is excited:
    # the shear branch (Ω ≈ 20 here against Ω ≈ 0.28) stays at zero amplitude,
    # which is precisely what makes the run a test — any of it that shows up in
    # the solution is discretisation error, and it is fast enough to be obvious.
    #-----------------------------------------------------------------------------
    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]
            x  = mesh.coords[1, ip]

            q.qn[ip, 1] = 0.0                             # ρA v
            q.qn[ip, 2] = 0.0                             # ρI ω
            q.qn[ip, 3] = (m.W*m.k - m.Φ)*cos(m.k*x)      # γ
            q.qn[ip, 4] = -m.Φ*m.k*sin(m.k*x)             # χ

            # No background state to subtract: the reference configuration of
            # this problem is the undeformed beam at rest.
            q.qe[ip, 1] = 0.0
            q.qe[ip, 2] = 0.0
            q.qe[ip, 3] = 0.0
            q.qe[ip, 4] = 0.0
        end
    end

    p = timoshenko_properties()
    println("   beam: L = ", p.L, ", h = ", p.h, ", L/h = ", p.L/p.h)
    println("   wave speeds: shear √(κₛG/ρ) = ", p.cs, ", extensional √(E/ρ) = ", p.ce)
    println("   mode n = ", m.n, ": Ω = ", m.Ω, ", period T = ", 2π/m.Ω)
    println(" Initialize fields for the 1D Timoshenko beam (simply supported) ... DONE ")

    return q
end
