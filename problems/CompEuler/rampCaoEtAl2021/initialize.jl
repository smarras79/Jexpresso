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
    # Uniform free stream everywhere at t = 0, i.e. an impulsive start.
    # The wall is not initialised to its own temperature: the no-slip
    # isothermal condition of user_bc.jl imposes it from the first stage,
    # and letting the boundary layer grow out of the free stream is the
    # cheapest way to a converged 2D field.  The first few hundred steps
    # are the strongest transient of the whole run -- see the note on
    # :Delta_t in user_inputs.jl.
    #
    for ip = 1:mesh.npoin
        q.qn[ip,1]   = ρ∞
        q.qn[ip,2]   = ρ∞*u∞
        q.qn[ip,3]   = ρ∞*v∞
        q.qn[ip,4]   = ρE∞
        q.qn[ip,end] = p∞

        # Reference state = the free stream.  This case runs in TOTAL()
        # mode, but qe is what the perturbation output and the PERT()
        # branches of the user_* routines read, and it is also the state
        # the DynSGS norms are taken relative to, so it is filled with a
        # meaningful state rather than zeros.
        q.qe[ip,1]   = ρ∞
        q.qe[ip,2]   = ρ∞*u∞
        q.qe[ip,3]   = ρ∞*v∞
        q.qe[ip,4]   = ρE∞
        q.qe[ip,end] = p∞
    end

    if rank == 0
        println(" Initialize fields for 2D CompEuler (rampCaoEtAl2021) ... DONE ")
    end

    return q
end
