#---------------------------------------------------------------------------------
# 3D emergence of a twisted convection-zone flux rope into a coronal point dipole,
# the parameter study of Leake, Linton & Antiochos, "The role of reconnection in
# the onset of solar eruptions", ApJ 934:10 (2022). See README.md.
#
#   julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022
#   mpiexec -n 128 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022
#
# :lk22_sim picks the simulation: the rope is rotated by θ = n π/4 about z (n = 0..7,
# the paper's simulations 0-7). n = 4 (default) is the eruption the paper describes
# in detail (its Figs. 3-4, 12); n = 7, 0, 1 are stable; 2 and 6 the intermediate
# cases. :lk22_dipole => false gives the paper's "No dipole" run (use it with n = 0).
#
# LK22_REF picks the mesh (LK22.geo) and the time step:
#   1:  16 x 16 x 24 elements  (smoke test)
#   2:  32 x 32 x 48 elements  (1/4 of the paper's resolution at :nop => 4)
#   3:  64 x 64 x 96 elements  (1/2; generate it first: gmsh -3 LK22.geo -setnumber ref 3 -o LK22_64x64x96.msh)
#   4: 128 x 128 x 192 elements (the paper's 512 x 512 x 768 cells; ref 4, LK22_128x128x192.msh)
#---------------------------------------------------------------------------------

# Positivity floors as the Runge-Kutta stage limiter (as problems/MHD/fluxEmergenceLeake2014):
# ρ ≥ f ρ_e(z) and p ≥ f p_e(z) of the reference atmosphere, p raised through E at fixed
# velocity and field. Not conservative where it fires; the hits are counted over all ranks and
# reported by rank 0 each time the total passes a power of ten.
if !@isdefined(lk22_floor_frac)
    const lk22_floor_frac = Ref{Float64}(0.2)
    const lk22_floor_hits = Ref{Int}(0)
    const lk22_floor_next = Ref{Int}(1)
end

# CFL monitor (the time step is fixed, the wave speeds are not): every 100 stages the largest
# (|v| + c_f) Δt/Δh over the nodes, with c_f = √((γp + |B|²)/ρ) and Δh the smallest LGL gap of
# the elements around the node, is reported by rank 0 when it has grown by 20 %, as a
# WARNING above 0.5 (the step is set for ≈ 0.25 initially: README).
if !@isdefined(lk22_cfl_dh)
    const lk22_cfl_dh    = Ref{Vector{Float64}}(Float64[])
    const lk22_cfl_count = Ref{Int}(0)
    const lk22_cfl_shown = Ref{Float64}(0.0)
end

function lk22_node_spacing(mesh)
    ngl = Int(mesh.ngl)
    h   = fill(Inf, Int(mesh.npoin))
    for ie = 1:Int(mesh.nelem)
        he = Inf
        for m = 1:ngl-1, (a, b) in ((mesh.connijk[ie,m,1,1], mesh.connijk[ie,m+1,1,1]),
                                    (mesh.connijk[ie,1,m,1], mesh.connijk[ie,1,m+1,1]),
                                    (mesh.connijk[ie,1,1,m], mesh.connijk[ie,1,1,m+1]))
            he = min(he, sqrt((mesh.coords[1,b] - mesh.coords[1,a])^2 + (mesh.coords[2,b] - mesh.coords[2,a])^2 +
                              (mesh.coords[3,b] - mesh.coords[3,a])^2))
        end
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[ie,i,j,k]
            h[ip] = min(h[ip], he)
        end
    end
    return h
end

function lk22_positivity_limiter!(u, integrator, p, t)
    npoin = p.mesh.npoin
    qe    = p.qp.qe
    frac  = lk22_floor_frac[]
    γm1   = γ_mhd - 1.0
    nhit  = 0
    lk22_cfl_count[] += 1
    docfl = lk22_cfl_count[] % 100 == 1
    if docfl && length(lk22_cfl_dh[]) != npoin
        lk22_cfl_dh[] = lk22_node_spacing(p.mesh)
    end
    dh    = lk22_cfl_dh[]
    Δt    = Float64(p.Δt)
    cfl = 0.0; vmx = 0.0; cmx = 0.0
    @inbounds for ip = 1:npoin
        iρ = ip; iρu = npoin + ip; iρv = 2npoin + ip; iE = 3npoin + ip; iρw = 4npoin + ip
        iBx = 5npoin + ip; iBy = 6npoin + ip; iBz = 7npoin + ip; iψ = 8npoin + ip
        ρmin = frac*qe[ip,1]
        pmin = frac*qe[ip,end]
        ρ = u[iρ]
        if ρ < ρmin
            u[iρ] = ρmin; ρ = ρmin; nhit += 1
        end
        ke = 0.5*(u[iρu]^2 + u[iρv]^2 + u[iρw]^2)/ρ
        me = 0.5*(u[iBx]^2 + u[iBy]^2 + u[iBz]^2) + 0.5*u[iψ]^2
        if γm1*(u[iE] - ke - me) < pmin
            u[iE] = pmin/γm1 + ke + me; nhit += 1
        end
        if docfl
            B2 = u[iBx]^2 + u[iBy]^2 + u[iBz]^2
            v  = sqrt(u[iρu]^2 + u[iρv]^2 + u[iρw]^2)/ρ
            cf = sqrt(max(γ_mhd*γm1*(u[iE] - ke - me) + B2, 0.0)/ρ)
            c  = (v + cf)*Δt/dh[ip]
            if c > cfl
                cfl = c; vmx = v; cmx = cf
            end
        end
    end
    if docfl
        comm = MPI.COMM_WORLD
        g = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(cfl, MPI.MAX, comm) : cfl
        if g > 1.2*lk22_cfl_shown[]
            if MPI.Comm_rank(comm) == 0       # the |v|, c_f printed are those of rank 0's largest
                print(string(g > 0.5 ? " # WARNING CFL: " : " # CFL: ", "max (|v| + c_f) Δt/Δh = ", round(g, sigdigits = 3),
                             " at t = ", round(t, digits = 3), " (rank 0: |v| = ", round(vmx, sigdigits = 3), ", c_f = ",
                             round(cmx, sigdigits = 3), ")", g > 0.5 ? "; reduce :Δt" : "", "\n"))
                flush(stdout)
            end
            lk22_cfl_shown[] = g
        end
    end
    # One report for all ranks, from rank 0, when the total passes 1, 10, 100, ...: every rank
    # calls the stage limiter at every Runge-Kutta stage, so the sum is a safe collective. t is
    # the stage time t_n + c_i Δt (in t₀), not an output time.
    comm = MPI.COMM_WORLD
    nall = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(nhit, +, comm) : nhit
    lk22_floor_hits[] += nall
    if nall > 0 && lk22_floor_hits[] >= lk22_floor_next[]
        if MPI.Comm_rank(comm) == 0
            print(string(" # positivity floors: ", lk22_floor_hits[], " node hits so far on all ranks (",
                         nall, " in this stage) at t = ", round(t, digits = 3), "\n"))
            flush(stdout)
        end
        while lk22_floor_next[] <= lk22_floor_hits[]
            lk22_floor_next[] *= 10
        end
    end
    return nothing
end

# Explicit dissipation (Sec. 2.1: η = 1.458e1 Ω m and ν = 3.350e3 kg m⁻¹ s⁻¹, "both equal
# to 0.01" in normalized units), per element centroid, added to the DynSGS coefficient of the
# momentum (ν, kinematic: README) and field (η, uniform) slots.
lk22_dissipation(x, y, z) = (0.01, 0.01)

function user_inputs()

    LK22_REF = 2
    mesh     = ("LK22_16x16x24.msh", "LK22_32x32x48.msh", "LK22_64x64x96.msh", "LK22_128x128x192.msh")[LK22_REF]
    # Δt for (|v| + c_f) Δt/Δh ≈ 0.25 at |v| + c_f ≈ 80 on the smallest LGL gap Δh of each mesh
    # (0.173 of the 8, 4, 2, 1 L₀ tall bottom elements): the fast speed just above the transition
    # region is 45 at t = 0 and reconnection outflows reach tens of v₀ (README). The CFL monitor
    # reports the actual value.
    Δt       = (4.0e-3, 2.0e-3, 1.0e-3, 5.0e-4)[LK22_REF]

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(stage_limiter! = lk22_positivity_limiter!,
                                                      williamson_condition = false),
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => 400.0,      # Fig. 8: simulation 4 erupts at 150-200 t₀, confined by 400
        # every 10 t₀, plus the times of the paper's Figs. 3-4 (142, 184); README for denser output
        :diagnostics_at_times => sort(vcat(collect(0.0:10.0:400.0), [142.0, 184.0])),
        :conservation_every   => 10,
        :conservation_slots   => [1, 4],
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,       # gravity, ψ damping, velocity damping
        :lwell_balanced       => true,       # subtract the discrete RHS of qe (atmosphere + dipole)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # The paper's problem (README):
        #---------------------------------------------------------------------------
        :lk22_sim             => 4,          # θ = n π/4: 4 erupts (Figs. 3-4); 7, 0, 1 stable; 2, 6 intermediate
        :lk22_dipole          => true,       # false: the paper's "No dipole" run (with :lk22_sim => 0)
        :lk22_rope            => true,       # false: the background alone, which must stay at rest
        :lk22_d               => 15000.0,    # dipole strength (moment along −x): 19.5 G at the photosphere
        :lk22_zd              => -100.0,     # dipole depth (not printed in the paper; README)
        :lk22_cz_gradient     => 1.0,        # Eq. 23 as printed (exactly adiabatic with m_f = 1)
        :lk22_mf              => 1.0,        # p = ρT/m_f: 1 reproduces the paper's Fig. 1 (README); 1.25 = literal text
        :lk22_glm_alpha       => 1.0,        # ψ damping rate α_p c_h/Δh
        :lk22_damp_x_frac     => 0.75,       # side damping from |x|, |y| = 3L_x/8 (Eq. 21)
        :lk22_damp_z_frac     => 0.875,      # top damping from z = 7 max(z)/8 (Eq. 22)
        :mhd_dissipation      => lk22_dissipation,
        #---------------------------------------------------------------------------
        # DynSGS in conserved form on the departure from the reference state, as
        # problems/MHD/fluxEmergenceLeake2014.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model           => DSGS_MHD(),
        :dsgs_sensor          => "legacy",
        :dsgs_CR              => 1.0,
        :dsgs_Cmax            => 0.5,
        :dsgs_Cmin            => 0.0,
        :dsgs_gamma           => 5.0/3.0,
        :dsgs_Prt             => 1.0,
        :dsgs_norms           => "element",
        :dsgs_local_rel       => 1.0,
        :dsgs_conserved       => true,
        :dsgs_fast_floors     => true,
        :dsgs_thermal_residual => true,
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh: LK22.geo (the paper's stretching, element boundaries every 32, 16, 8, 4 cells)
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/fluxEmergenceLeake2022/" * mesh,
        :lfilter              => false,
        :outformat            => "vtk",
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :output_dir           => "./output/",
        :loutput_pert         => false,
        :linitial_refine      => false,
        :ladapt               => false,
    ) #Dict

    return inputs
end
