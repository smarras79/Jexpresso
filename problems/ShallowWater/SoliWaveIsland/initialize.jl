function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    """
    SoliWaveIsland: solitary wave run-up/run-down on a circular conical island.

    Sec. 5.5 of Marras, Kopera, Constantinescu, Suckale, Giraldo (2018).

    Initial state:
      Synolakis (1987) solitary surface profile centred at xc_wave, propagating
      rightward over a still water layer of depth h0:

          η(x, 0)  = A · sech²( γ · (x - xc_wave) )
          γ        = sqrt(3 A / (4 h0³))
          c        = sqrt( g (h0 + η) )
          u(x, 0)  = c η / (h0 + η)
          v(x, 0)  = 0

      Bathymetry (cone):
          Hb(x,y) = 0.93 · max(0, 1 - r/rc),   r = √((x-xc_cone)² + (y-yc_cone)²)

      State vector q = [H, Hu, Hv] where H is the local water depth above the
      bathymetry, so that

          H(x, y, 0) = max(H_dry, h0 + η - Hb)
    """

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D NLSWE solitary-wave-on-island (Marras et al. 2018, Sec. 5.5) ..."
    end

    qvars    = ["H", "Hu", "Hv"]
    qoutvars = ["H", "Hu", "Hv"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat,
                 inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if (inputs[:backend] == CPU())

        # --- physical constants ---------------------------------------------------
        g        = 9.81

        # --- solitary wave (Synolakis 1987) --------------------------------------
        A        = 0.064   # amplitude [m]
        h0       = 0.32    # still water depth [m]
        xc_wave  = 2.5     # initial wave centre [m]
        γ        = sqrt(3.0 * A / (4.0 * h0^3))

        # --- conical island ------------------------------------------------------
        xc_cone  = 12.5
        yc_cone  = 0.0
        rc_cone  = 3.6
        hc_cone  = 0.93

        # --- wet/dry safety floor ------------------------------------------------
        H_dry    = 1.0e-4

        if rank == 0
            @info "  h0   = $(h0) m,   A = $(A) m,  γ = $(γ) 1/m"
            @info "  cone: centre=($xc_cone, $yc_cone),  rc=$rc_cone,  hmax=$hc_cone"
        end

        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            y = mesh.y[ip]

            # surface displacement at t=0 (sech² profile)
            arg   = γ * (x - xc_wave)
            sech_ = 1.0 / cosh(arg)
            η     = A * sech_ * sech_

            # bathymetry (cone) and water depth
            dx = x - xc_cone
            dy = y - yc_cone
            r  = sqrt(dx*dx + dy*dy)
            Hb = r < rc_cone ? hc_cone * (1.0 - r / rc_cone) : 0.0

            H_water = h0 + η - Hb
            H_water = max(H_water, H_dry)

            # rightward-propagating solitary wave: u = c η / (h0 + η)
            if r < rc_cone
                # over the island, suppress incident velocity to avoid initial impulse
                u_x = 0.0
            else
                c   = sqrt(g * (h0 + η))
                u_x = c * η / (h0 + η)
            end

            q.qn[ip, 1] = H_water
            q.qn[ip, 2] = H_water * u_x
            q.qn[ip, 3] = 0.0

            # reference / "exact" state used by the kernel: still water at rest
            H_eq = max(h0 - Hb, H_dry)
            q.qe[ip, 1] = H_eq
            q.qe[ip, 2] = 0.0
            q.qe[ip, 3] = 0.0
        end
    end

    if rank == 0
        @info " Initialize fields for 2D NLSWE solitary-wave-on-island ... DONE"
    end

    return q
end
