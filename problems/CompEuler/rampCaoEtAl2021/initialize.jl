#---------------------------------------------------------------------------------
# 2D hypersonic compression-ramp flow with laminar separation.
#
#   S. Cao, J. Hao, I. Klioutchnikov, H. Olivier, C.-Y. Wen,
#   "Unsteady effects in a hypersonic compression ramp flow with laminar
#    separation", J. Fluid Mech. 912, A3 (2021), doi:10.1017/jfm.2020.1093
#
# The paper's DNS is three-dimensional and periodic in the span; this deck
# is the x-y version of it, i.e. the TWO-DIMENSIONAL flow that the paper
# itself computes twice -- as the starting field of the 3D runs (Section
# 2.3: "the three-dimensional simulations are initialised by duplicating
# the two-dimensional converged solution in the spanwise direction") and as
# the base flow of the global stability analysis of Section 3.2.  Nothing
# in the 2D problem is an approximation of the 3D one: the geometry is
# two-dimensional, and the three-dimensionality of the paper is an
# INSTABILITY of exactly this flow, not an ingredient of it.
#
# What this case reproduces (2D, steady once converged):
#   - the laminar flat-plate boundary layer, delta = 1.38 mm at separation
#   - the separation bubble, 0.59 < x/L < 1.26 (Section 3.1)
#   - the separation shock, the reattachment shock and their interaction
#   - the surface pressure and heat flux of figure 2(a) / figure 3(c)
#
# What it cannot reproduce, by construction: the streamwise heat-flux
# streaks, the spanwise modulation and the low-frequency unsteadiness of
# Sections 3-5.  Those are the global 3D instability of this base flow.
#
# GEOMETRY (Section 2.2, model of Roghelia et al. 2017b)
#   sharp leading edge, flat plate L = 100 mm, ramp 15 deg also 100 mm,
#   1 mm of free stream ahead of the leading edge.  See ramp15.geo.
#
# FLOW CONDITIONS (Table 1, shock tunnel TH2, RWTH Aachen)
#   M_inf = 7.7,  T_inf = 125 K,  p_inf = 760 Pa,  u_inf = 1726 m/s,
#   Re_inf = 4.2e6 1/m  (Re_inf,L = 4.2e5),  isothermal wall T_w = 293 K,
#   gamma = 1.4, Pr = 0.71, Sutherland's law for mu.
#
# The total enthalpy, 1.7 MJ/kg, is low enough for the perfect gas
# assumption the paper makes (Section 2.2), so no high-temperature
# chemistry is needed here either.
#---------------------------------------------------------------------------------

#
# Free-stream state, shared by initialize() and by the supersonic-inflow /
# free-stream boundary conditions in user_bc.jl.
#
# Returns (rho, u, v, p, T, rhoE) in SI units.
#
# Table 1 quotes M = 7.7, T = 125 K, p = 760 Pa and u = 1726 m/s.  The
# velocity is rebuilt here from M_inf and the code's own PhysConst rather
# than hardcoded, so the stream is exactly Mach 7.7 for the gas the solver
# integrates; what comes out is 1725.6 m/s, i.e. the paper's value.
#
# CONSISTENCY OF THE REYNOLDS NUMBER.  The paper gives BOTH Re_inf and the
# conditions it is built from, so the two must agree, and they pin down the
# viscosity law:
#
#   rho_inf = p/(R T)                            = 0.021185 kg/m^3
#   mu(125 K) by Sutherland (1.716e-5, 273.15, 110.4)
#                                                = 8.656e-6 Pa s
#   Re_inf,L = rho_inf u_inf L / mu(T_inf)       = 4.22e5
#
# against the paper's 4.2e5 -- so the standard air Sutherland constants,
# which is what user_inputs.jl sets, are the ones the paper used.  Change
# :sutherland_muref and this case stops being the paper's case.
#
function ramp_freestream()

    PhysConst = PhysicalConst{Float64}()

    M∞ = 7.7                                  # free stream Mach number
    p∞ = 760.0                                # Pa
    T∞ = 125.0                                # K

    ρ∞ = p∞/(PhysConst.Rair*T∞)               # kg/m^3
    c∞ = sqrt(PhysConst.γ*PhysConst.Rair*T∞)  # m/s
    u∞ = M∞*c∞
    v∞ = 0.0

    ρE∞ = p∞/PhysConst.γm1 + 0.5*ρ∞*(u∞*u∞ + v∞*v∞)

    return ρ∞, u∞, v∞, p∞, T∞, ρE∞
end

#
# Isothermal wall temperature, Table 1.  Used by user_bc.jl.
#
ramp_Twall() = 293.0   # K


#---------------------------------------------------------------------------------
# STARTING FIELD: a compressible laminar boundary layer, not a uniform stream.
#
# WHY NOT A UNIFORM FREE STREAM.  problems/CompEuler/ffs_step starts from
# one, and copying that here does not work -- it is the single reason this
# case failed on its first run, at step 3.  ffs_step's walls are FREE SLIP,
# so a uniform stream already satisfies its boundary conditions at t = 0.
# This wall is NO SLIP and ISOTHERMAL, so a uniform stream violates it at
# every wall node: u steps 1726 -> 0 and T steps 125 -> 293 across the
# first LGL interval, 8e-6 m, along the whole plate and ramp.  That is a
# step discontinuity INSIDE a spectral element, and the discretisation has
# no way to carry it -- DynSGS least of all, since its residual sensor is
# held at zero for the first `:dsgs_hold_steps + 1` steps (rhs.jl,
# _dsgs_hold_steps) precisely so that it cannot misread a fresh initial
# condition.  The solution is unstabilised exactly when it is least
# defensible, and it dies on the step the sensor wakes up.
#
# So the starting field carries the boundary layer from the first instant:
#
#   velocity   Pohlhausen quartic u/ue = 2z - 2z^3 + z^4 in the Howarth
#              variable, i.e. the laminar shape with the right wall slope
#              and a smooth join to the free stream at z = 1
#   temperature Crocco-Busemann, T = Tw + (Taw - Tw)(u/ue) + (Te - Taw)(u/ue)^2,
#              which is exactly Tw at the wall and Te at the edge, with the
#              recovery bump (peak ~530 K here) in between
#   density    from the boundary-layer approximation dp/dn = 0, so p = p_inf
#              across the layer and rho = p_inf/(R T)
#
# The wall-normal coordinate is stretched by the Howarth integral
# int (T/Te) dz, which is what makes the cold dense sublayer thin and the
# hot outer part fat -- get this wrong and the profile has the right
# thickness but the wrong wall gradient, which is the one thing that
# matters for the heat flux.
#
# delta(s) is calibrated on the paper's own number, delta = 1.38 mm at the
# separation location s = 59 mm (Section 2.2), and grown as sqrt(s) from
# there.  It is floored at 1.4e-4 m -- three wall elements -- over the
# first 0.6 mm, because the sharp leading edge is a genuine singularity
# where the true delta is thinner than the grid and no starting field can
# resolve it.
#
# This is NOT the converged solution: it carries no shock, no separation,
# no pressure rise.  It is a smooth, boundary-condition-consistent field
# that the run develops FROM, and it is consistent to the last digit where
# it matters -- at the wall it gives rho = p_inf/(R*293) and hence
# rho*E = p_inf/(gamma-1), which is exactly what the isothermal no-slip
# condition in user_bc.jl imposes.  Section 2.5 / figure 2(b) compares the
# CONVERGED profiles against the similarity solution; this is the same
# family of profile used as a starting guess, not as an answer.
#---------------------------------------------------------------------------------

#
# Wall-normal profile, tabulated once: returns (zeta_table, y_over_delta).
#
function ramp_profile_table(Te, Tw, Taw; nz=400)
    z  = collect(range(0.0, 1.0, length=nz+1))
    Tt = similar(z)
    for i in eachindex(z)
        su    = 2z[i] - 2z[i]^3 + z[i]^4          # u/ue, Pohlhausen
        Tt[i] = Tw + (Taw - Tw)*su + (Te - Taw)*su*su
    end
    yy = similar(z); yy[1] = 0.0
    for i = 1:nz
        yy[i+1] = yy[i] + 0.5*(Tt[i] + Tt[i+1])/Te*(z[i+1] - z[i])   # Howarth
    end
    yy ./= yy[end]
    return z, yy
end

#
# (u/ue, T) at a wall-normal distance n inside a layer of thickness delta.
#
function ramp_profile_at(n, δ, z, yy, Te, Tw, Taw)
    (n >= δ || δ <= 0.0) && return 1.0, Te
    target = n/δ
    lo, hi = 1, length(yy)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        yy[mid] < target ? (lo = mid) : (hi = mid)
    end
    w  = (yy[hi] - yy[lo]) > 0 ? (target - yy[lo])/(yy[hi] - yy[lo]) : 0.0
    ζ  = z[lo] + w*(z[hi] - z[lo])
    su = 2ζ - 2ζ^3 + ζ^4
    return su, Tw + (Taw - Tw)*su + (Te - Taw)*su*su
end


function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        println(" Initialize fields for 2D CompEuler (rampCaoEtAl2021: Mach-7.7 15-deg compression ramp) ... ")
    end

    #---------------------------------------------------------------------------------
    # Solution variables.
    #
    # Slot 4 is the TOTAL ENERGY, not rho*theta: theta is an entropy
    # variable, conserved across a contact but not across a shock, so the
    # Euler-theta system carries the wrong speed for the separation and
    # reattachment shocks no matter how it is stabilized.  :energy_equation
    # => "energy" in user_inputs.jl selects the matching kernel branches.
    #
    # NOTICE: the length of qvars defines neqs; qoutvars can hold at most
    # neqs+1 entries, which is why the output stops at T.
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρE"]
    qoutvars = ["ρ", "u", "v", "p", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    ρ∞, u∞, v∞, p∞, T∞, ρE∞ = ramp_freestream()

    if rank == 0
        PhysConst = PhysicalConst{Float64}()
        M∞  = u∞/sqrt(PhysConst.γ*p∞/ρ∞)
        μ∞  = 1.716e-5*(T∞/273.15)^1.5*(273.15 + 110.4)/(T∞ + 110.4)
        ReL = ρ∞*u∞*0.1/μ∞
        @printf("    free stream: rho = %.6f kg/m3, u = %.1f m/s, p = %.1f Pa, T = %.1f K, M = %.3f\n",
                ρ∞, u∞, p∞, T∞, M∞)
        @printf("    mu(T_inf) = %.4e Pa.s -> Re_inf,L = %.3e  (paper: 4.2e5)\n", μ∞, ReL)
        @printf("    isothermal wall T_w = %.1f K,  T_w/T_inf = %.2f\n", ramp_Twall(), ramp_Twall()/T∞)
    end

    #
    # Geometry of the wall, mirroring ramp15.geo: sharp leading edge at the
    # origin, plate to x = L, then a ramp at alpha.  For a node we need the
    # arclength s from the leading edge and the wall-normal distance n.
    #
    PhysConst = PhysicalConst{Float64}()
    L, α = 0.1, deg2rad(15.0)
    cα, sα = cos(α), sin(α)

    Tw  = ramp_Twall()
    r   = sqrt(0.71)                                  # laminar recovery factor
    Taw = T∞*(1.0 + r*0.5*PhysConst.γm1*7.7^2)        # ~1374 K
    zt, yyt = ramp_profile_table(T∞, Tw, Taw)

    δref, sref = 1.38e-3, 0.059                       # Section 2.2
    #
    # delta(s) IS SMOOTH AND THE FLOOR IS RESOLVABLE.  Two defects met in one
    # cell, and a 1-rank run of rampCaoEtAl2021_M7 named it exactly: the
    # global first negative pressure was at (x, y) = (7.103784e-4,
    # 1.496259e-4), the mid-LGL node of streamwise element 2, at the top of
    # the first wall element, y/delta = 0.988.
    #
    # (1) max(delta_floor, delta_ref*sqrt(s/sref)) has a SLOPE DISCONTINUITY
    #     where the branches cross, at s = sref*(delta_floor/delta_ref)^2 =
    #     6.0723e-4 m -- 28% of the way along that same element. A
    #     C0-but-not-C1 field inside a spectral element is a Gibbs generator,
    #     and at M = 7.7 the internal energy is 5.7% of the total, so the
    #     ringing reaches p amplified 17.5x. sqrt(delta_floor^2 +
    #     delta_ref^2*s/sref) has the same two asymptotes with no kink.
    #
    # (2) 1.4e-4 was ONE element height (0.06/401 = 1.496e-4) on
    #     ramp15_uniform.msh, so the whole starting boundary layer was five
    #     LGL nodes thick there against twelve elements at x = 100 mm -- the
    #     leading edge was 12x less resolved than the rest of the plate. The
    #     old comment claimed "~3 wall elements" while setting one.
    #
    # 6.0e-4 is four element heights, and it dominates over the first
    # sref*(6.0e-4/1.38e-3)^2 = 1.115e-2 m = 11 mm of plate. Separation is at
    # 59 mm, so the flow this case exists to compute is untouched; what is
    # smeared is the sharp-leading-edge singularity, which no starting field
    # on any grid resolves anyway.
    #
    δfloor     = 6.0e-4                               # 4 element heights (0.06/401)

    nbl = 0
    for ip = 1:mesh.npoin
        # @view(mesh.coords[1,:])/@view(mesh.coords[2,:]) are deprecated; the node coordinates live in
        # mesh.coords[dim, ip] (3 x npoin) on the current kernel.
        x, y = mesh.coords[1,ip], mesh.coords[2,ip]

        if x <= 1.0e-12
            s_wall, n_wall = -1.0, y                  # the leading-edge column
        elseif x <= L
            s_wall, n_wall = x, y                     # flat plate
        else
            ξ = x - L                                 # ramp frame
            s_wall = L + ξ*cα + y*sα
            n_wall = -ξ*sα + y*cα
        end

        if s_wall <= 0.0 || n_wall <= 0.0
            u, v, T = u∞, v∞, T∞                      # free stream
        else
            δ  = sqrt(δfloor*δfloor + δref*δref*(s_wall/sref))   # smooth: no kink
            su, T = ramp_profile_at(n_wall, δ, zt, yyt, T∞, Tw, Taw)
            n_wall < δ && (nbl += 1)
            #-------------------------------------------------------------------
            # DIRECTION: wall-tangent at the wall, HORIZONTAL at the edge of
            # the layer.  The turn belongs to the boundary layer, not to the
            # free stream.
            #
            # THE BUG THIS REPLACES.  ramp_profile_at returns su = 1 for every
            # node OUTSIDE the layer, so
            #
            #     u, v = su*u∞*cα, su*u∞*sα        for every node with x > L
            #
            # turned the UNDISTURBED FREE STREAM by 15 degrees over the whole
            # block above the ramp -- up to the top boundary 60 mm away, where
            # nothing turns the flow before the shock exists.  At t = 0:
            #
            #   * along the entire vertical line x = 0.1, from the wall to
            #     y = 60 mm, v jumped 0 -> u∞ sin15 = 446.4 m/s; and
            #   * along the entire "top" boundary above the ramp, user_bc.jl
            #     prescribes the free stream, v = 0, while the node one LGL
            #     interval below it carried 446.4 m/s.
            #
            # A 446 m/s shear across 2.58e-5 m held open by a Dirichlet
            # condition -- the same illegal starting field this file's header
            # warns about for the no-slip wall, at the TOP of the domain.
            # MEASURED on rampCaoEtAl2021_M7: the global first positivity
            # repair was at (x, y) = (0.10036588, 0.06007220) on RHS call 180,
            # i.e. 0.37 mm past the ramp corner and exactly one LGL interval
            # below the top boundary, on step 36 -- when the fastest signal
            # had travelled 0.07 mm and the wall was 60 mm away.  Nothing
            # propagated there; this line put it there.
            #
            # theta(n) = alpha*(1 - su) is wall-tangent where su = 0 and
            # horizontal where su = 1.  Continuous in n, reduces to the old
            # plate behaviour at alpha = 0, and outside the layer it is the
            # EXACT free stream, so it agrees with the inflow and top
            # Dirichlet conditions instead of fighting them.
            #-------------------------------------------------------------------
            θ    = (x <= L) ? 0.0 : α*(1.0 - su)
            u, v = su*u∞*cos(θ), su*u∞*sin(θ)
        end

        # Boundary-layer approximation: p is constant across the layer.
        p = p∞
        ρ = p/(PhysConst.Rair*T)
        ρE = p/PhysConst.γm1 + 0.5*ρ*(u*u + v*v)

        q.qn[ip,1]   = ρ
        q.qn[ip,2]   = ρ*u
        q.qn[ip,3]   = ρ*v
        q.qn[ip,4]   = ρE
        q.qn[ip,end] = p

        # Reference state = the FREE STREAM, not the starting field.  qe is
        # what the DynSGS norms measure the departure from and what the
        # perturbation output subtracts, and both want the undisturbed
        # stream as the datum.
        q.qe[ip,1]   = ρ∞
        q.qe[ip,2]   = ρ∞*u∞
        q.qe[ip,3]   = ρ∞*v∞
        q.qe[ip,4]   = ρE∞
        q.qe[ip,end] = p∞
    end

    #
    # GLOBAL, not rank-local.  The rank-local form printed "laminar BL on 0 of
    # 56333 nodes" on 32 ranks, which reads as "this starting field has no
    # boundary layer"; in fact rank 0's partition holds no near-wall nodes.
    # Both sums count shared interface nodes once per owning rank, so npo_g
    # slightly exceeds the true node count -- but nbl_g is summed the same
    # way, so the ratio is right and that is what the line is for.
    #
    nbl_g = MPI.Allreduce(nbl,        MPI.SUM, comm)
    npo_g = MPI.Allreduce(mesh.npoin, MPI.SUM, comm)
    if rank == 0
        @printf("    starting field: laminar BL on %d of %d nodes (global; %d of %d on rank 0), T_aw = %.0f K\n",
                nbl_g, npo_g, nbl, mesh.npoin, Taw)
    end

    if rank == 0
        println(" Initialize fields for 2D CompEuler (rampCaoEtAl2021) ... DONE ")
    end

    return q
end
