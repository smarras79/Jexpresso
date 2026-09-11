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
#   JEXPRESSO_SV_TEND   final time                              (default 1.0)
#   JEXPRESSO_SV_CMIN   the DynSGS background floor :dsgs_Cmin  (default 0)
#   JEXPRESSO_SV_VISC   "dsgs" (default) or "none" for the plain Galerkin run,
#                       the two panels of the paper's Fig. 1
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
_sv_visc() = lowercase(strip(get(ENV, "JEXPRESSO_SV_VISC", "dsgs")))

function _sv_dt()
    d = tryparse(Float64, get(ENV, "JEXPRESSO_SV_DT", ""))
    d === nothing || return d
    return SV_DT_REF*64.0/max(1, _sv_nelx()*_sv_nop())
end

_sv_mesh() = string("./problems/MHD/smoothVortex/vortex_", _sv_nelx(), "x", _sv_nelx(), ".msh")

function user_inputs()
    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
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
        :dsgs_sensor      => "residual",
        :dsgs_CR          => 1.0,
        :dsgs_Cmax        => 0.5,
        :dsgs_Cmin        => _sv_cmin(),   # no background floor: it would be an
                                           # order-independent O(h) error on a
                                           # smooth solution (see the README)
        :dsgs_gamma       => 5.0/3.0,
        :dsgs_Prt         => 1.0,
        :dsgs_conserved   => true,
        :dsgs_norms       => "rank",
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
