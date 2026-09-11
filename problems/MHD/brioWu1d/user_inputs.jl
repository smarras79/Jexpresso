#---------------------------------------------------------------------------------
# Brio-Wu MHD shock tube, 1D ideal MHD (Brio & Wu 1988), as set up in Dao &
# Nazarov, J. Sci. Comput. 92:77 (2022), Sec. 5.2: domain (0, 1), γ = 2,
# t ∈ [0, 0.1] (see below), left/right states in initialize.jl. The MHD counterpart of
# problems/CompEuler/sod1d: stabilized by DynSGS (kernel/physics/SGS.jl,
# the 1D DSGS_MHD kernel) in its conserved form, integrated with SSPRK53.
#
#   julia --project=. src/Jexpresso.jl MHD brioWu1d
#
# 150 elements at N = 4 are 600 LGL points, the "600 DOFs" of the paper's
# Fig. 2(a); :nelx => 300 gives its 1200-DOF case. Final time 0.1 on (0, 1)
# = Brio & Wu's 0.2 on (−1, 1), the state the paper's Fig. 2 shows.
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Polynomial order of this run. `JEXPRESSO_BW_NOP` overrides the deck so the
# multi-order comparison of user_plot.jl can be produced without editing any
# file, and so that every order is run at the SAME number of degrees of
# freedom (~600, the paper's Fig. 2(a) resolution): unless the element count
# is given explicitly with `JEXPRESSO_BW_NELX`, it is chosen as 600/nop.
#
#     tools/brio_wu_order_scan.sh        # the whole sweep
#
# or by hand,
#
#     for D in 150 300 600 1200; do
#       for N in 4 5 6 7; do
#         JEXPRESSO_BW_DOFS=$D JEXPRESSO_BW_NOP=$N \
#             julia --project=. src/Jexpresso.jl MHD brioWu1d
#       done
#     done
#
# The time step follows the resolution, Δt = 5e-5 * (600/DOFs), so the Courant
# and the DynSGS parabolic numbers are the same at every point of the sweep
# and the comparison is not contaminated by a changing time error. At 600 DOFs
# the smallest LGL spacing is 1.15e-3 at nop 4 (150 elements) and 7.5e-4 at
# nop 7 (86 elements), so the Courant number against the fast speed of the
# right state runs from 0.16 to 0.25 and the parabolic number stays below 0.25.
#
# Overrides, all optional:
#   JEXPRESSO_BW_NOP    polynomial order                       (default 4)
#   JEXPRESSO_BW_DOFS   degrees of freedom, sets the element count (default 600)
#   JEXPRESSO_BW_NELX   element count, overrides JEXPRESSO_BW_DOFS
#   JEXPRESSO_BW_DT     time step, overrides the rule above
#   JEXPRESSO_BW_CMIN   the DynSGS background floor :dsgs_Cmin (default 0.06)
#---------------------------------------------------------------------------------
const BW_DOFS_DEFAULT = 600
const BW_DT_AT_DEFAULT = 5.0e-5

_bw_nop()  = something(tryparse(Int,     get(ENV, "JEXPRESSO_BW_NOP",  "")), 4)
_bw_dofs() = something(tryparse(Int,     get(ENV, "JEXPRESSO_BW_DOFS", "")), BW_DOFS_DEFAULT)
_bw_cmin() = something(tryparse(Float64, get(ENV, "JEXPRESSO_BW_CMIN", "")), 0.06)

# The coefficient per NODE (Dao & Nazarov's own form, their eq. 4.10) instead
# of per element. At a fixed number of degrees of freedom the element form
# gives the higher orders fewer, larger elements and therefore a coarser
# staircase of nu, which is not a fair comparison across orders; the nodal
# form has no such bias.
_bw_nodal() = lowercase(strip(get(ENV, "JEXPRESSO_BW_NODAL", "false"))) in ("1", "true", "yes", "on")

# Element length scale of the DynSGS coefficient. The kernels use
# Δ = Δ_K/(k+1) (`ngl`); Dao & Nazarov's eq. 4.10 uses h_K/k, the polynomial
# DEGREE. The difference is a factor (k+1)/k, and it is not innocent when
# orders are compared at a FIXED number of degrees of freedom: with
# n_elx = DOFs/k,
#
#     Δ_K/(k+1) = (L/DOFs)·k/(k+1)      grows with k  (0.8 at k=4, 0.875 at k=7)
#     Δ_K/k     =  L/DOFs               the same for every order
#
# so the (k+1) convention gives the higher orders a systematically larger
# viscosity, while the paper's is order-independent. Δ enters the cap and the
# floor linearly and the residual viscosity quadratically, so the convention
# can be switched here without touching the kernels, by scaling the
# coefficients: C_max, C_min by (k+1)/k and C_R by ((k+1)/k)².
#
# THE ORDER COMPARISON OF THIS CASE USES THE PAPER'S Δ_K/k BY DEFAULT: it is
# the only one of the two that is order-independent at fixed DOFs, and the
# comparison exists to reproduce their Fig. 2. Nothing outside this deck
# changes — `src/kernel/physics/SGS.jl` still computes Δ = Δ_K/ngl for every
# case, and every other deck (orszagTangBormanis2024, sod1d, thetaTracers, the
# forward-facing step, the rising bubble) is untouched by this switch.
#
#     JEXPRESSO_BW_HSCALE=nop     the paper's Δ_K/k          (the default here)
#     JEXPRESSO_BW_HSCALE=ngl     Δ_K/(k+1), the kernels' own, unscaled
const BW_HSCALE_DEFAULT = "nop"

_bw_hscale() = lowercase(strip(get(ENV, "JEXPRESSO_BW_HSCALE", BW_HSCALE_DEFAULT)))

function _bw_hfac()
    _bw_hscale() == "nop" || return 1.0
    N = _bw_nop()
    return (N + 1)/N
end

function _bw_nelx()
    n = tryparse(Int, get(ENV, "JEXPRESSO_BW_NELX", ""))
    n === nothing || return n
    return max(1, round(Int, _bw_dofs()/_bw_nop()))
end

# Δt ∝ 1/DOFs: the Courant and parabolic numbers are then the same at every
# resolution of a convergence sweep.
function _bw_dt()
    d = tryparse(Float64, get(ENV, "JEXPRESSO_BW_DT", ""))
    d === nothing || return d
    return BW_DT_AT_DEFAULT*BW_DOFS_DEFAULT/max(1, _bw_nelx()*_bw_nop())
end

function user_inputs()
    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        # Fastest signal: the fast magnetosonic speed of the right state,
        # √((γp + |B|²)/ρ) = 3.75, against the smallest LGL spacing
        # 0.146·(1/150) ≈ 1e-3: Δt = 5e-5 is a Courant number of 0.19.
        :Δt                   => _bw_dt(),   # 5e-5 at 600 DOFs, scaled by 1/DOFs
        :tinit                => 0.0,
        # t = 0.1 on (0, 1) is Brio & Wu's t = 0.2 on (−1, 1), the state of
        # the paper's Fig. 2 (its "t̂ = 0.2"): see README.md.
        :tend                 => 0.1,
        :diagnostics_at_times => (0:0.025:0.1),
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        :lperiodic_1d         => false,
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => _bw_nop(),   # JEXPRESSO_BW_NOP overrides; see the top of this file
        #---------------------------------------------------------------------------
        # DynSGS-MHD, conserved form: one residual-based kinematic ν per
        # element (max over the equations of the normalized residual, Dao &
        # Nazarov eq. 4.8; :dsgs_CR, :dsgs_Cmax are their C_R, C_max), capped at
        # C_max·Δ·(|u| + c_f), applied as ∇·(ν∇q) to every conserved variable
        # (user_primitives.jl). :μ are the per-slot multipliers for
        # (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz); :dsgs_gamma must equal γ_mhd = 2
        # of user_flux.jl.
        #---------------------------------------------------------------------------
        :lvisc            => true,
        :μ                => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model       => DSGS_MHD(),
        :dsgs_sensor      => "residual",  # element-wise strong residual (DSGS.md §1.2); "legacy" = the pre-Sep-2026 sensor
        :dsgs_CR          => 1.0*_bw_hfac()^2,
        :dsgs_Cmax        => 0.5*_bw_hfac(),
        # Background floor C_min (not in the paper), 6 % of the first-order
        # viscosity C_max·h·(|u|+c_f). (0.03 was enough before the residual's
        # time derivative was made stage-consistent, DSGS.md §4.4: the old
        # stencil added its own dissipation on the moving waves.)
        # The residual viscosity is C_R h² R: it scales with the amplitude
        # of what it sees, so the element-scale ripples that the slowly
        # moving compound wave radiates into the plateau behind it (±0.5 %
        # in ρ, one wiggle per element, present at P3 and P4, with the
        # element and the nodal coefficient alike) are never damped by it.
        # The floor damps them at a rate ν(π/h)² ≈ 300/unit time while
        # diffusing a resolved profile by √(2νt) ≈ 0.004 over the run; 0.01
        # leaves a trace of them, 0.03 none (measured). Dao & Nazarov's P3
        # elements with exact quadrature do not show these ripples; the
        # collocated LGL flux of this code is the remaining difference.
        :dsgs_Cmin        => _bw_cmin()*_bw_hfac(),  # JEXPRESSO_BW_CMIN overrides
        :dsgs_gamma       => 2.0,
        :dsgs_Prt         => 1.0,
        :dsgs_conserved   => true,
        # Nodal form of the kernel (Dao & Nazarov's own): ν at every node from
        # the assembled residual, so the coefficient is a continuous field and
        # the diffusive flux has no jump at element interfaces; the element
        # form (one ν per element) left one wiggle per element in the plateau
        # behind the compound wave. :dsgs_Cl is their C_l = 0.4 (eq. 4.7).
        # With the C_min floor both forms give the same clean profile; the
        # element form (the default) is kept here.
#       :ldsgs_nodal      => true,    # false (the default) = one ν per element
        :dsgs_Cl          => 0.4,
        :ldsgs_nodal      => _bw_nodal(),  # JEXPRESSO_BW_NODAL=1 for the paper's nodal form
        :dsgs_norms       => "rank",    # residual normalized by the spread over this rank's elements (default; "domain": the whole tube, MPI-reduced)
        :energy_equation  => "energy",
        :lkep              => false,
        :entropy_variables => false,
        :lfilter           => false,
        #---------------------------------------------------------------------------
        # Mesh: uniform 1D grid built from xmin/xmax/nelx
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => 0.0,
        :xmax                 => 1.0,
        :nelx                 => _bw_nelx(),  # 600/nop, so every order runs at ~600 DOFs
        #---------------------------------------------------------------------------
        # Output. Default: the density figure of the paper's Fig. 2
        # (user_plot.jl), density-it<n>.png, with the reference solution of
        # reference_hll.dat and its three zoom boxes at the final time.
        # :plot_user => false gives instead the format of CompEuler/sod1d: one
        # figure fields-it<n>.png per output time with a panel per output
        # variable (ρ, u, v, p, By; the reference dashed at the final time)
        # and a last panel with the DynSGS coefficient per element
        # (:plot_matrix => false writes those panels as separate files).
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :plot_user            => true,
        :loverwrite_output    => true,
        :output_dir           => "./output",
    )
    return inputs
end
