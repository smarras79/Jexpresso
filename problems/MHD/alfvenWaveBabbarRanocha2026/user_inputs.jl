#---------------------------------------------------------------------------------
# Smooth Alfvén wave: the GLM-MHD convergence test of Babbar & Ranocha,
# J. Comput. Phys. 563 (2026) 115060, §5.5.1 (Tóth, J. Comput. Phys. 161
# (2000) 605). Exact solution in initialize.jl; [0, √2]², doubly periodic,
# γ = 5/3, run to t = 2 (two periods, so the exact solution at the end is the
# initial condition).
#
# The paper's Fig. 13 is the L² error at t = 2 against the number of degrees of
# freedom for N = 1, 2, 3 on 16², 32², 64², 128² elements. `user_plot.jl` here
# measures the same thing at the final time (errors/nop<N>_nelx<M>_<visc>.dat)
# and draws every stored error as that figure. Jexpresso is run at N >= 3, so
# the sweep is N = 3 (the paper's highest order) to 7.
#
#     tools/alfven_wave_scan.sh        # the sweep: N = 3 ... 7, DSGS and Galerkin, and its figure
#     tools/alfven_wave_plot.jl        # the figure and the table again, from the stored errors
#     tools/alfven_wave_mesh.sh        # the meshes (4² ... 128² are in the repo)
#
# or one run by hand,
#
#     JEXPRESSO_AW_NOP=3 JEXPRESSO_AW_NELX=32 \
#         julia --project=. src/Jexpresso.jl MHD alfvenWaveBabbarRanocha2026
#
# Overrides, all optional:
#   JEXPRESSO_AW_NOP    polynomial order N                      (default 3)
#   JEXPRESSO_AW_NELX   elements per side; picks the mesh file  (default 16)
#   JEXPRESSO_AW_DT     time step, overrides the rule below
#   JEXPRESSO_AW_TEND   final time                              (default 2.0)
#   JEXPRESSO_AW_VISC   "dsgs" (default) or "none" for the plain Galerkin run
#   JEXPRESSO_AW_SENSOR :dsgs_sensor, "residual" (default) | "legacy"
#   JEXPRESSO_AW_CMIN   :dsgs_Cmin, the DynSGS background floor (default 0)
#   JEXPRESSO_AW_FREEZE 1 or 0 forces :dsgs_freeze_stage (the coefficient once
#                       per step, or at every stage); unset, mod_inputs.jl's
#                       default applies: on with the residual sensor
#   JEXPRESSO_AW_SOLVER ck54, vern7, vern9, dp8, ssprk54 — the time integrator;
#                       unset, the rule below (ck54 up to N = 5, vern7 above)
#
# The time step follows the resolution: a Courant number of 0.1 on the mean
# node spacing √2/(nelx N) against the fastest wave of the initial condition
# (|v| + c_f <= 0.1 + sqrt(γp + |B|²) = 1.185), rounded down so that a whole
# number of steps reaches t_end. On the smallest LGL spacing that is 0.12 for
# N = 3, 0.14 for N = 4 and 0.22 for N = 7.
#
# THE TIME INTEGRATOR follows N, so that its error stays below the spatial one
# at that Courant number. CarpenterKennedy2N54 is 4th order: its stability
# polynomial is 1 + z + z²/2 + z³/6 + z⁴/24 + z⁵/200, so it shifts the phase of
# the wave (ω = 2π) by (ωΔt)⁵/300 per step, and by t = 2 that is an error of
# about 9 Δt⁴ in the state (the wave's L² norm is 0.14). For N = 7 on 16² that
# is 1.2e-11 against a spatial error of 5.5e-12 (measured: 1.32e-11 with it,
# 5.53e-12 with Vern7 or Vern9); for N = 6 on 32² it would be a third of the
# spatial error (1.4e-12 against 3.75e-12). Vern7 (7th order) is used from N = 6 on; it costs 10 right-hand sides
# per step against CK54's 5. The density error, the paper's measure, does not
# see this error at all: ρ = 1 whatever the phase of the wave.
#---------------------------------------------------------------------------------
const AW_NELX_DEFAULT = 16
const AW_NOP_DEFAULT  = 3
const AW_COURANT      = 0.1
const AW_LAMBDA       = 0.1 + sqrt(5.0/3.0*0.1 + 1.01)   # |v| + c_f of the initial condition
const AW_LBOX         = sqrt(2.0)                           # 1/cos(π/4) = 1/sin(π/4)

_aw_nop()  = something(tryparse(Int,     get(ENV, "JEXPRESSO_AW_NOP",  "")), AW_NOP_DEFAULT)
_aw_nelx() = something(tryparse(Int,     get(ENV, "JEXPRESSO_AW_NELX", "")), AW_NELX_DEFAULT)
_aw_tend() = something(tryparse(Float64, get(ENV, "JEXPRESSO_AW_TEND", "")), 2.0)
_aw_cmin() = something(tryparse(Float64, get(ENV, "JEXPRESSO_AW_CMIN", "")), 0.0)
_aw_visc() = lowercase(strip(get(ENV, "JEXPRESSO_AW_VISC", "dsgs")))
_aw_sensor() = lowercase(strip(get(ENV, "JEXPRESSO_AW_SENSOR", "residual")))
_aw_freeze() = (v = strip(get(ENV, "JEXPRESSO_AW_FREEZE", ""));
                isempty(v) ? nothing : v in ("1", "true", "yes"))

function _aw_dt()
    d = tryparse(Float64, get(ENV, "JEXPRESSO_AW_DT", ""))
    d === nothing || return d
    dt   = AW_COURANT*AW_LBOX/(max(1, _aw_nelx()*_aw_nop())*AW_LAMBDA)
    tend = _aw_tend()
    return tend/ceil(tend/dt)
end

# the integrator of the rule above for order N
_aw_default_solver(N) = N <= 5 ? "ck54" : "vern7"

function _aw_solver_named(name)
    name == "vern7"   && return Vern7()
    name == "vern9"   && return Vern9()
    name == "dp8"     && return DP8()
    name == "ssprk54" && return SSPRK54()
    name == "ck54"    || @warn "alfvenWaveBabbarRanocha2026: unknown JEXPRESSO_AW_SOLVER=$(name); using CarpenterKennedy2N54"
    return CarpenterKennedy2N54()
end

function _aw_solver()
    name = lowercase(strip(get(ENV, "JEXPRESSO_AW_SOLVER", "")))
    return _aw_solver_named(isempty(name) ? _aw_default_solver(_aw_nop()) : name)
end

_aw_mesh() = string("./problems/MHD/alfvenWaveBabbarRanocha2026/alfven_",
                    _aw_nelx(), "x", _aw_nelx(), ".msh")

function user_inputs()
    inputs = Dict(
        :ode_solver           => _aw_solver(),
        :Δt                   => _aw_dt(),
        :tinit                => 0.0,
        :tend                 => _aw_tend(),
        :diagnostics_at_times => (0.0, _aw_tend()),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,   # GLM ψ-damping source (Dedner mixed cleaning)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => _aw_nop(),
        #---------------------------------------------------------------------------
        # Stabilization: DynSGS as set up for the smooth vortex, the other
        # accuracy test of this directory (residual sensor, no background
        # floor, domain-wide normalization, physical form). On a smooth
        # solution the residual is small and so is ν; JEXPRESSO_AW_VISC=none
        # gives the plain Galerkin run to compare against.
        #---------------------------------------------------------------------------
        :lvisc            => (_aw_visc() != "none"),
        :μ                => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model       => DSGS_MHD(),
        :dsgs_sensor      => _aw_sensor(),
        :dsgs_CR          => 1.0,
        :dsgs_Cmax        => 0.5,
        :dsgs_Cmin        => _aw_cmin(),   # 0: a floor would add an O(h) error
        :dsgs_gamma       => 5.0/3.0,
        :dsgs_Prt         => 1.0,
        :dsgs_conserved      => false,     # physical form, as user_primitives.jl
        :dsgs_nazarov_energy => true,
        :dsgs_norms       => "domain",     # independent of the rank count
        :lrichardson      => false,
        :energy_equation  => "energy",
        :lkep              => false,
        :entropy_variables => false,
        :lfilter           => false,
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => _aw_mesh(),
        #---------------------------------------------------------------------------
        :outformat           => "png",
        :plot_user           => true,
        :plot_matrix         => false,
        :plot_vars           => ["ρ", "u", "w", "By", "Bz", "p"],
        :loverwrite_output   => true,
        :lwrite_initial      => false,
        :output_dir          => "./output",
        :loutput_pert        => false,
        :linitial_refine     => false,
        :ladapt              => false,
    ) #Dict
    f = _aw_freeze()
    f === nothing || (inputs[:dsgs_freeze_stage] = f)
    return inputs
end
