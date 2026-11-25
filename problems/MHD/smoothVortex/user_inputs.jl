#---------------------------------------------------------------------------------
# Smooth MHD vortex — the accuracy test of Dao & Nazarov,
# J. Sci. Comput. 92:77 (2022), §5.1.1 (Balsara, ApJS 151:149 (2004)).
#
# A steady vortex on a uniform advection over the doubly periodic
# [-L/2, L/2]² (L = 20 by default, the paper's [-10,10]²), so the exact
# solution is the initial condition translated by v₀t and the error can be
# measured exactly. This is the case their Fig. 1 convergence history is built
# on, and `user_plot.jl` here reproduces that figure: the L¹, L² and L∞ error
# of the velocity against #DOFs on log-log axes, one line per polynomial order
# and stabilization, with slope guides.
#
#     tools/smooth_vortex_mpi_scan.sh      # the sweep, DynSGS and Galerkin
#                                          # (README.md: the options)
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
#   JEXPRESSO_SV_VISC   "dsgs" (default) or "none" for the plain Galerkin run,
#                       the two panels of the paper's Fig. 1
#   JEXPRESSO_SV_L      the box width L; picks vortex_L<L>_<n>x<n>.msh (default 20)
#
# Every DynSGS setting (sensor, norms, C_1, C_max, C_min, rel, cutoff, hold,
# BDF order, freeze) is a literal entry of the Dict below; edit it there. The
# environment variables that used to override them are refused (below).
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
# The DynSGS settings are literal entries of the Dict in user_inputs(). These
# environment variables used to override them and are no longer read: a run
# that still sets one stops, instead of silently running the deck's value.
const SV_DSGS_ENV_REMOVED = ("JEXPRESSO_SV_SENSOR", "JEXPRESSO_SV_NORMS", "JEXPRESSO_SV_C1",
                             "JEXPRESSO_SV_CR", "JEXPRESSO_SV_CMAX", "JEXPRESSO_SV_CMIN",
                             "JEXPRESSO_SV_REL", "JEXPRESSO_SV_CUTOFF", "JEXPRESSO_SV_HOLD",
                             "JEXPRESSO_SV_BDF", "JEXPRESSO_SV_FREEZE")
_sv_visc() = lowercase(strip(get(ENV, "JEXPRESSO_SV_VISC", "dsgs")))

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
# The default is the box of the paper, [-10,10]^2, i.e. L = 20.
_sv_lbox() = something(tryparse(Float64, get(ENV, "JEXPRESSO_SV_L", "")), 20.0)
# The box is ALWAYS in the mesh name (vortex_L20_32x32.msh): a mesh whose name
# does not say which box it is cannot be told apart from one that is a
# different box, and reading the wrong one costs a whole sweep.
_sv_ltag() = (L = _sv_lbox();
                 string("L", L == round(L) ? string(Int(round(L))) : string(L), "_"))
_sv_mesh() = string("./problems/MHD/smoothVortex/vortex_", _sv_ltag(),
                       _sv_nelx(), "x", _sv_nelx(), ".msh")

function user_inputs()
    for v in SV_DSGS_ENV_REMOVED
        haskey(ENV, v) && error(" problems/MHD/smoothVortex: ", v, " is no longer read; set the DynSGS settings in user_inputs.jl.")
    end
    inputs = Dict(
        :ode_solver           => _sv_solver(),
        :Δt                   => _sv_dt(),
        :tinit                => 0.0,
        :tend                 => _sv_tend(),
        :diagnostics_at_times => (0.0, _sv_tend()),
        # Mass and energy, the discrete totals Σ_K Σ_ij ω_iω_j J q of ρ and ρE,
        # at every step, as for the rotor; tools/plot_conservation.py draws
        # |Q(t) − Q(0)|/|Q(0)|. One file per run, named as its error record
        # (user_plot.jl, _sv_tag), since the runs of a sweep share the output
        # directory: conservation_nop<N>_nelx<M>_<tag>.dat.
        :conservation_every   => 1,
        :conservation_slots   => [1, 4],
        :conservation_file    => inp -> string("conservation_nop", _sv_nop(), "_nelx", _sv_nelx(),
                                               "_", _sv_tag(inp), ".dat"),
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
        # "residual": each element's own strong residual against the BDF
        # history (an inter-element mismatch reaches the sensor; the assembled
        # rate of the former "legacy" sensor cancels it).
        :dsgs_sensor      => "residual",
        :dsgs_rel         => 1.0,   # floor = the variable's own physical scale (1e-3: the pre-Sep-2026 rule)
        :dsgs_C1          => 1.0,
        :dsgs_Cmax        => 0.5,
        :dsgs_Cmin        => 0.0,          # no background floor: it would be an
                                           # order-independent O(h) error on a
                                           # smooth solution (see the README)
        :dsgs_gamma       => 5.0/3.0,
        :dsgs_Prt         => 1.0,
        # PHYSICAL form, the one user_primitives.jl is written for (ρ, u, v,
        # T = p/ρ, w, B, ψ): ν on ρ, ρ̄ν in the deviatoric stress, κ = ρ̄ν/Pr
        # on ∇T (Dao & Nazarov 2022 §4.4, as in orszagTangBormanis2024), η = ν
        # on B, and the τ·u and η B·∇B work terms in the energy. With
        # :dsgs_conserved => true these primitives got one kinematic ν each
        # (ν∇u, ν∇T, no stress, no work terms) — neither form.
        :dsgs_conserved      => false,
        :dsgs_nazarov_energy => true,
        # THE NORMALIZATION SCOPE CHANGES THE ANSWER, and by a lot. "rank"
        # normalizes each rank's residual by ITS OWN subdomain's spread, so ν
        # — and therefore the solution — depends on how the domain was cut.
        # Measured on this case (nop 6, 32² elements, Δt = 3.3333e-4, t = 1,
        # absolute velocity L¹): 3.935e-7 on 2 ranks against 7.155e-6 on 4,
        # an 18x difference from nothing but the partition, and a floor that
        # no mesh refinement can go below. "domain" costs two small Allreduces
        # per RHS call (a few doubles) and is what Dao & Nazarov's eq. 4.8
        # means by the norm over Ω.
        :dsgs_norms       => "domain",
        # The coefficient once per step instead of once per stage: it is what
        # keeps the residual, and with it nu, from stalling at O(dt) and
        # capping a high-order accuracy test.
        :dsgs_freeze_stage => true,
        # Smoothness cutoff on the normalized residual (SGS.jl, _dsgs_cut): ν
        # takes max(0, ratio − cutoff). With the coefficient frozen per step the
        # sensor still reads the BDF2 truncation error of a moving smooth
        # solution, (Δt²/3)∂ₜ³q, which alone gives ν ∝ Δ²Δt² ∝ h⁴ and caps the
        # measured order near 4 from N = 5 on (setup_notes/dsgsCutoff.pdf). A
        # shock's ratio is O(10²), so 1e-3 removes the floor and leaves a shock
        # as it was. Records: "dsgs" at 1e-3, "dsgs_nocut" at 0, "dsgs_cut"
        # otherwise (user_plot.jl, SV_DSGS_CUTOFF: keep that equal to this).
        :dsgs_cutoff       => 1.0e-3,
        # Steps the coefficient is held at zero at the start of a run. 2 is the
        # minimum the BDF2 history needs; more probes whether the smooth-flow
        # excess is a STARTUP dose (records "dsgs_hold").
        :dsgs_hold_steps   => 2,
        # Order k of the residual sensor's backward difference: with BDF-k the
        # floor above is −Δtᵏ/(k+1)·∂ₜᵏ⁺¹q (records "_bdf<k>" for k ≠ 2).
        :dsgs_bdf_order    => 2,
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
