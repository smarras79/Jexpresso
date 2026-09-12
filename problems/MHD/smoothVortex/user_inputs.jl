#---------------------------------------------------------------------------------
# Smooth MHD vortex — the accuracy test of Dao & Nazarov,
# J. Sci. Comput. 92:77 (2022), §5.1.1 (Balsara, ApJS 151:149 (2004)).
#
# A steady vortex on a uniform advection over the doubly periodic [-5,5]²,
# so the exact solution is the initial condition translated by v₀t and the
# error can be measured exactly. This is the case their Fig. 1 convergence
# history is built on, and `user_plot.jl` here reproduces that figure:
# the L¹, L² and L∞ error of the velocity against 1/sqrt(#DOFs) on log-log
# axes, one line per polynomial order, with slope guides.
#
#     tools/smooth_vortex_mesh.sh          # the meshes, once
#     tools/smooth_vortex_scan.sh          # the sweep: 4 orders x 4 meshes
#
# or by hand,
#
#     JEXPRESSO_SV_NOP=3 JEXPRESSO_SV_NELX=16 \
#         julia --project=. src/Jexpresso.jl MHD smoothVortex
#
# Overrides, all optional:
#   JEXPRESSO_SV_NOP    polynomial order                        (default 4)
#   JEXPRESSO_SV_NELX   elements per side; picks the mesh file  (default 16)
#   JEXPRESSO_SV_DT     time step, overrides the rule below
#   JEXPRESSO_SV_SOLVER ck54 (default), vern9, vern7, dp8, ssprk54, tsit5 —
#                       the time integrator; see the note by _sv_solver on why
#                       a 4th-order one caps the measured rate at 4
#   JEXPRESSO_SV_TEND   final time                              (default 1.0)
#   JEXPRESSO_SV_CMIN   the DynSGS background floor :dsgs_Cmin  (default 0)
#   JEXPRESSO_SV_CR     :dsgs_CR                                 (default 1)
#   JEXPRESSO_SV_REL    :dsgs_rel, the normalization floor        (default 1)
#   JEXPRESSO_SV_NORMS  :dsgs_norms, "domain" (default) | "rank" | "element";
#                       "rank" makes the answer depend on the rank count
#   JEXPRESSO_SV_CMAX   :dsgs_Cmax                               (default 0.5)
#   JEXPRESSO_SV_VISC   "dsgs" (default) or "none" for the plain Galerkin run,
#                       the two panels of the paper's Fig. 1
#   JEXPRESSO_SV_SENSOR "residual" (default, the element's own residual) or
#                       "legacy" (the assembled rate)
#
# The time step follows the resolution, Δt ∝ 1/(nelx·nop), so the Courant
# number is the same at every point of a sweep and the comparison is not
# contaminated by a changing time error. The reference value below is a
# Courant number of about 0.2 against the fastest wave of this initial
# condition (|v| + c_f ≈ 3.1) on the smallest LGL spacing.
#---------------------------------------------------------------------------------
const SV_NELX_DEFAULT = 16
const SV_NOP_DEFAULT  = 4
const SV_DT_REF       = 2.0e-3    # at nelx*nop = 64

_sv_nop()  = something(tryparse(Int,     get(ENV, "JEXPRESSO_SV_NOP",  "")), SV_NOP_DEFAULT)
_sv_nelx() = something(tryparse(Int,     get(ENV, "JEXPRESSO_SV_NELX", "")), SV_NELX_DEFAULT)
_sv_tend() = something(tryparse(Float64, get(ENV, "JEXPRESSO_SV_TEND", "")), 1.0)
_sv_cmin() = something(tryparse(Float64, get(ENV, "JEXPRESSO_SV_CMIN", "")), 0.0)
# C_R and C_max as well, so that a run can keep the SENSOR (and its
# diagnostics) while applying no viscosity at all: with both zero the
# solution is the plain Galerkin one and the printed residual is the
# residual of a clean solution, which is how one checks that the sensor
# itself converges under refinement.
_sv_cr()   = something(tryparse(Float64, get(ENV, "JEXPRESSO_SV_CR",   "")), 1.0)
# :dsgs_rel — how far below its own physical scale a variable's spread may fall
# before that scale normalizes its residual. 1 is the fix; 1e-3 is what the
# kernels used before it, and reproduces the first-order convergence this case
# was built to expose.
_sv_rel()  = something(tryparse(Float64, get(ENV, "JEXPRESSO_SV_REL",  "")), 1.0)
# "domain" (the default here) | "rank" | "element" — see :dsgs_norms below.
_sv_norms() = lowercase(strip(get(ENV, "JEXPRESSO_SV_NORMS", "domain")))
_sv_cmax() = something(tryparse(Float64, get(ENV, "JEXPRESSO_SV_CMAX", "")), 0.5)
_sv_visc() = lowercase(strip(get(ENV, "JEXPRESSO_SV_VISC", "dsgs")))
# "residual" (the default) measures the ELEMENT's own strong residual against
# the BDF history; "legacy" measures the ASSEMBLED rate instead. Which one is
# used decides whether an inter-element mismatch can reach the sensor, so the
# convergence test needs to be able to run both.
_sv_sensor() = lowercase(strip(get(ENV, "JEXPRESSO_SV_SENSOR", "residual")))

function _sv_dt()
    d = tryparse(Float64, get(ENV, "JEXPRESSO_SV_DT", ""))
    d === nothing || return d
    return SV_DT_REF*64.0/max(1, _sv_nelx()*_sv_nop())
end

#---------------------------------------------------------------------------------
# THE TIME ERROR IS PART OF A CONVERGENCE SWEEP, and on this case it is what
# limits the measured rate. CarpenterKennedy2N54 is FOURTH order, and the rule
# above takes Δt ∝ h, so the total error is
#
#       C_s h^{N+1} + C_t Δt^4 ∝ h^{N+1} + h^4,
#
# and no order above 3 can show its own rate: p saturates at 4 however fine
# the mesh. Measured on the plain Galerkin runs (nop 4, 4/8/16 elements per
# side): 5.47, then 4.12, then 3.64 as the second term takes over.
#
# Two ways out, both switchable here:
#   JEXPRESSO_SV_DT=<fixed>      one Δt for the whole sweep (tools/
#                                smooth_vortex_scan.sh sets it from the finest
#                                run of the sweep), so the time error is a
#                                CONSTANT and stops polluting the slope until
#                                it dominates the finest point;
#   JEXPRESSO_SV_SOLVER=vern9    an eighth/ninth-order integrator, so the time
#                                error is below the spatial one at any Δt this
#                                case can run.
#
# `vern9` (Vern9, 9th order) and `dp8` (DP8, 8th) are the two worth having;
# `ck54` is the default low-storage CarpenterKennedy2N54, `ssprk54` the
# strong-stability-preserving 5-stage 4th-order one the AdvDiff cases use.
#---------------------------------------------------------------------------------
function _sv_solver()
    name = lowercase(strip(get(ENV, "JEXPRESSO_SV_SOLVER", "ck54")))
    name == "vern9"   && return Vern9()
    name == "vern7"   && return Vern7()
    name == "dp8"     && return DP8()
    name == "ssprk54" && return SSPRK54()
    name == "tsit5"   && return Tsit5()
    name == "ck54"    || @warn "smoothVortex: unknown JEXPRESSO_SV_SOLVER=$(name); using CarpenterKennedy2N54"
    return CarpenterKennedy2N54()
end

# The BOX WIDTH. The vortex is a Gaussian, so on [-L/2, L/2]² the exact
# solution is not periodic: the velocity perturbation at the middle of an
# edge is (L/2)exp((1 - (L/2)²)/2)/2π with opposite sign on opposite edges,
# so the initial condition jumps across the periodic seam by twice that —
# 4.9e-06 on the L = 10 box of Balsara and of Dao & Nazarov. That jump is a
# discontinuity in the DATA: it floors any accuracy study at ~1e-5, whatever
# the order, and a 6th-order element reaches the floor sooner than a 4th.
# L = 15 puts the floor at 1e-12 and L = 20 at machine zero; generate those
# meshes with SV_L=20 tools/smooth_vortex_mesh.sh.
_sv_lbox() = something(tryparse(Float64, get(ENV, "JEXPRESSO_SV_L", "")), 10.0)
_sv_ltag() = (L = _sv_lbox(); L == 10.0 ? "" :
                 string("L", L == round(L) ? string(Int(round(L))) : string(L), "_"))
_sv_mesh() = string("./problems/MHD/smoothVortex/vortex_", _sv_ltag(),
                       _sv_nelx(), "x", _sv_nelx(), ".msh")

function user_inputs()
    inputs = Dict(
        :ode_solver           => _sv_solver(),
        :Δt                   => _sv_dt(),
        :tinit                => 0.0,
        :tend                 => _sv_tend(),
        :diagnostics_at_times => (0.0, _sv_tend()),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,   # GLM ψ-damping source (Dedner mixed cleaning)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => _sv_nop(),
        #---------------------------------------------------------------------------
        # Stabilization. The point of the test is that the residual viscosity
        # must NOT destroy the accuracy of a high-order solution on a smooth
        # problem, so the comparison of the paper's Fig. 1 is DynSGS against
        # the plain Galerkin run (JEXPRESSO_SV_VISC=none).
        #---------------------------------------------------------------------------
        :lvisc            => (_sv_visc() != "none"),
        :μ                => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model       => DSGS_MHD(),
        :dsgs_sensor      => _sv_sensor(),
        :dsgs_rel         => _sv_rel(),
        :dsgs_CR          => _sv_cr(),
        :dsgs_Cmax        => _sv_cmax(),
        :dsgs_Cmin        => _sv_cmin(),   # no background floor: it would be an
                                           # order-independent O(h) error on a
                                           # smooth solution (see the README)
        :dsgs_gamma       => 5.0/3.0,
        :dsgs_Prt         => 1.0,
        :dsgs_conserved   => true,
        # THE NORMALIZATION SCOPE CHANGES THE ANSWER, and by a lot. "rank"
        # normalizes each rank's residual by ITS OWN subdomain's spread, so ν
        # — and therefore the solution — depends on how the domain was cut.
        # Measured on this case (nop 6, 32² elements, Δt = 3.3333e-4, t = 1,
        # absolute velocity L¹): 3.935e-7 on 2 ranks against 7.155e-6 on 4,
        # an 18x difference from nothing but the partition, and a floor that
        # no mesh refinement can go below. "domain" costs two small Allreduces
        # per RHS call (a few doubles) and is what Dao & Nazarov's eq. 4.8
        # means by the norm over Ω.
        :dsgs_norms       => _sv_norms(),
        :lrichardson      => false,
        :energy_equation  => "energy",
        :lkep              => false,
        :entropy_variables => false,
        :lfilter           => false,
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => _sv_mesh(),
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
