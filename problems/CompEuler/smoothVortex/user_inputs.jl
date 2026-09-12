#---------------------------------------------------------------------------------
# Isentropic (Shu) vortex — the classical smooth accuracy test of the 2D
# compressible Euler equations, and the HYDRODYNAMIC CONTROL for
# problems/MHD/smoothVortex.
#
# Same box, same meshes, same error machinery and the same figures as the MHD
# case; no magnetic field, so no ∇·B and no GLM cleaning. The difference
# between the two convergence histories is the bias the divergence constraint
# introduces, which is what this case exists to measure.
#
#     tools/smooth_vortex_mesh.sh                          # the meshes, once
#     SV_CASE=CompEuler/smoothVortex tools/smooth_vortex_mpi_scan.sh
#
# or by hand,
#
#     JEXPRESSO_EV_NOP=4 JEXPRESSO_EV_NELX=16 \
#         julia --project=. src/Jexpresso.jl CompEuler smoothVortex
#
# Overrides, all optional:
#   JEXPRESSO_EV_NOP     polynomial order                        (default 4)
#   JEXPRESSO_EV_NELX    elements per side; picks the mesh file  (default 16)
#   JEXPRESSO_EV_DT      time step, overrides the rule below
#   JEXPRESSO_EV_TEND    final time                              (default 1.0)
#   JEXPRESSO_EV_BETA    vortex strength                         (default 5)
#                        5 is the classical test; 1 matches the amplitude of
#                        the MHD vortex (κ = μ = 1) for a like-for-like
#                        comparison with problems/MHD/smoothVortex
#   JEXPRESSO_EV_VISC    "dsgs" (default) or "none" for plain Galerkin
#   JEXPRESSO_EV_SOLVER  ck54 (default), vern9, vern7, dp8, ssprk54, tsit5
#   JEXPRESSO_EV_CMIN / _CR / _CMAX / _REL / _NORMS   the DynSGS coefficients
#
# THE MESHES ARE THE MHD CASE'S. They are the same [-5,5]² doubly periodic
# quad grids, and duplicating them here would only mean two copies to
# regenerate; the per-case SEM cache lives next to this deck, so nothing is
# shared but the .msh file itself.
#---------------------------------------------------------------------------------
const EV_NELX_DEFAULT = 16
const EV_NOP_DEFAULT  = 4
const EV_DT_REF       = 2.0e-3    # at nelx*nop = 64

_ev_nop()   = something(tryparse(Int,     get(ENV, "JEXPRESSO_EV_NOP",   "")), EV_NOP_DEFAULT)
_ev_nelx()  = something(tryparse(Int,     get(ENV, "JEXPRESSO_EV_NELX",  "")), EV_NELX_DEFAULT)
_ev_tend()  = something(tryparse(Float64, get(ENV, "JEXPRESSO_EV_TEND",  "")), 1.0)
_ev_beta()  = something(tryparse(Float64, get(ENV, "JEXPRESSO_EV_BETA",  "")), EV_BETA_DEFAULT)
_ev_cmin()  = something(tryparse(Float64, get(ENV, "JEXPRESSO_EV_CMIN",  "")), 0.0)
_ev_cr()    = something(tryparse(Float64, get(ENV, "JEXPRESSO_EV_CR",    "")), 1.0)
_ev_cmax()  = something(tryparse(Float64, get(ENV, "JEXPRESSO_EV_CMAX",  "")), 0.5)
_ev_rel()   = something(tryparse(Float64, get(ENV, "JEXPRESSO_EV_REL",   "")), 1.0)
_ev_visc()  = lowercase(strip(get(ENV, "JEXPRESSO_EV_VISC",  "dsgs")))
_ev_norms() = lowercase(strip(get(ENV, "JEXPRESSO_EV_NORMS", "domain")))

# The time integrator. As on the MHD case, a fourth-order one is not what
# limits these runs — see the note by _sv_solver in
# problems/MHD/smoothVortex/user_inputs.jl for the measurement.
function _ev_solver()
    name = lowercase(strip(get(ENV, "JEXPRESSO_EV_SOLVER", "ck54")))
    name == "vern9"   && return Vern9()
    name == "vern7"   && return Vern7()
    name == "dp8"     && return DP8()
    name == "ssprk54" && return SSPRK54()
    name == "tsit5"   && return Tsit5()
    name == "ck54"    || @warn "CompEuler/smoothVortex: unknown JEXPRESSO_EV_SOLVER=$(name); using CarpenterKennedy2N54"
    return CarpenterKennedy2N54()
end

# Δt ∝ 1/(nelx·nop) unless it is given: the Courant number is then the same at
# every point of a sweep. A convergence sweep should fix it instead
# (tools/smooth_vortex_mpi_scan.sh does), so that the time error is a constant
# rather than something that shrinks with h and contaminates the slope.
function _ev_dt()
    d = tryparse(Float64, get(ENV, "JEXPRESSO_EV_DT", ""))
    d === nothing || return d
    return EV_DT_REF*64.0/max(1, _ev_nelx()*_ev_nop())
end

_ev_mesh() = string("./problems/MHD/smoothVortex/vortex_", _ev_nelx(), "x", _ev_nelx(), ".msh")

function user_inputs()
    inputs = Dict(
        :ode_solver           => _ev_solver(),
        :Δt                   => _ev_dt(),
        :tinit                => 0.0,
        :tend                 => _ev_tend(),
        :diagnostics_at_times => (0.0, _ev_tend()),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => false,   # the vortex is an exact solution, unforced
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => _ev_nop(),
        #---------------------------------------------------------------------------
        # Stabilization. The point of the test is that the residual viscosity
        # must NOT destroy the accuracy of a high-order solution on a smooth
        # problem, so the comparison is DynSGS against the plain Galerkin run
        # (JEXPRESSO_EV_VISC=none).
        #---------------------------------------------------------------------------
        :lvisc            => (_ev_visc() != "none"),
        :μ                => [1.0, 1.0, 1.0, 1.0],
        :visc_model       => DSGS(),
        :dsgs_sensor      => "residual",
        :dsgs_CR          => _ev_cr(),
        :dsgs_Cmax        => _ev_cmax(),
        :dsgs_Cmin        => _ev_cmin(),   # no background floor: it would be an
                                           # order-independent O(h) error on a
                                           # smooth solution (see the README)
        :dsgs_Prt         => 0.7,
        :dsgs_rel         => _ev_rel(),
        :dsgs_norms       => _ev_norms(),
        :lrichardson      => false,
        :energy_equation  => "energy",
        :lkep             => false,
        :entropy_variables => false,
        :lfilter          => false,
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => _ev_mesh(),
        #---------------------------------------------------------------------------
        :outformat           => "png",
        :plot_user           => true,
        :plot_matrix         => false,
        :plot_vars           => ["ρ", "u", "v", "p"],
        :loverwrite_output   => true,
        :lwrite_initial      => false,
        :output_dir          => "./output",
        :loutput_pert        => false,
        :linitial_refine     => false,
        :ladapt              => false,
    ) #Dict
    return inputs
end
