function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the 2D plane-strain cantilever ............. ")

    #-----------------------------------------------------------------------------
    # U = (ρu, ρv, σxx, σyy, σxy, uₓ, u_y)ᵀ
    #
    # Rows 1-5 are the elastodynamic conservation laws; rows 6-7 are the
    # passive displacement recovery (zero flux, u̇ = v in the source) that turns
    # the output into a picture of the beam actually bending — warp the mesh by
    # (uₓ, u_y) in ParaView.
    #
    # qoutvars divide the momenta through by ρ so the VTK files carry the
    # particle velocity itself.
    #-----------------------------------------------------------------------------
    qvars    = ["ρu", "ρv", "σxx", "σyy", "σxy", "ux", "uy"]
    qoutvars = ["u",  "v",  "σxx", "σyy", "σxy", "ux", "uy"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    P  = elastic_properties()
    Ω1 = cantilever_omega1()

    #-----------------------------------------------------------------------------
    # Struck, not pulled: the beam starts UNDEFORMED and UNSTRESSED, with a
    # transverse velocity in the shape of the first bending mode,
    #
    #     u = 0,   v = V Φ₁(x),   σ = 0,   uₓ = u_y = 0
    #
    # V is chosen so the tip swings to ≈1.5% of L, which keeps the run well
    # inside linear elasticity: for simple harmonic motion the displacement
    # amplitude is V/Ω₁, so V = 0.015 Ω₁.
    #
    # WHY THIS AND NOT A RELEASED STATIC DEFLECTION. Every boundary condition of
    # this case is satisfied exactly by this state at t = 0: σ ≡ 0 makes the
    # three free surfaces traction-free, and Φ₁(0) = Φ₁'(0) = 0 makes the clamp
    # motionless. Releasing a statically loaded beam instead means switching off
    # a tip traction at t = 0, which is a step in time at that face and launches
    # a front the non-dissipative scheme has no way to absorb.
    #
    # Φ₁ is a mode of the reduced 1D beam model, not of the 2D elastic beam, so
    # the motion will carry a little of the neighbouring modes too. That is not
    # error — it is the difference between the two models, and it is the reason
    # to run this in 2D at all.
    #-----------------------------------------------------------------------------
    V = 0.015*Ω1

    for ip = 1:mesh.npoin
        x = mesh.x[ip]

        q.qn[ip, 1] = 0.0                       # ρu
        q.qn[ip, 2] = P.ρ*V*cantilever_mode1(x) # ρv
        q.qn[ip, 3] = 0.0                       # σxx
        q.qn[ip, 4] = 0.0                       # σyy
        q.qn[ip, 5] = 0.0                       # σxy
        q.qn[ip, 6] = 0.0                       # uₓ
        q.qn[ip, 7] = 0.0                       # u_y

        # No background state: the reference configuration is the undeformed
        # beam at rest.
        for ieq = 1:7
            q.qe[ip, ieq] = 0.0
        end
    end

    println("   material: ρ = ", P.ρ, ", E = ", P.E, ", ν = ", P.ν,
            "  →  λ = ", P.λ, ", μ = ", P.μ)
    println("   wave speeds: cp = ", P.cp, ", cs = ", P.cs)
    println("   beam: L = ", P.L, ", h = ", P.h, ", L/h = ", P.L/P.h,
            " (plane strain: bending modulus E/(1-ν²) = ", P.Ebend, ")")
    println("   mode 1 (Euler-Bernoulli estimate): Ω₁ = ", Ω1, ", period = ", 2π/Ω1)
    println("   tip velocity V = ", V, "  →  tip deflection ≈ ", V/Ω1)
    println(" Initialize fields for the 2D plane-strain cantilever ............. DONE ")

    return q
end
