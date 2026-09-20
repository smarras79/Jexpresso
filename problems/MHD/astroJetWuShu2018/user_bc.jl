#---------------------------------------------------------------------------------
# Boundary conditions — Wu & Shu (2018), Example 5.6:
#
#   "The fixed inflow beam condition is specified on the nozzle
#    {y = 0, |x| < 0.05}, and the others are outflow boundary conditions."
#
# The `tag` argument is the gmsh physical-curve name carried by the edge, so
# the four names below are exactly the groups declared in AJ.geo:
#
#   "bottom"  y = 0      the nozzle for |x| <= 0.05, outflow outside it
#   "right"   x = +0.5   outflow
#   "top"     y = 1.5    outflow
#   "left"    x = -0.5   outflow
#
# OUTFLOW. build_custom_bcs_dirichlet! (kernel/boundaryconditions/BCs.jl)
# pre-fills qbdy with a sentinel and copies back only the slots this routine
# overwrites, so simply NOT touching qbdy is the open/do-nothing condition: the
# interior solution convects out through the strong-form CG boundary untouched,
# and RHS is left alone there. That is what "outflow" means here. It is the
# same treatment as problems/CompEuler/ffs_step's supersonic outflow and it is
# only characteristically correct where the outgoing flow is supersonic; on the
# bottom boundary outside the nozzle, where the cocoon eventually pushes gas
# back down at subsonic speed, it is weakly reflective. That is tolerable here
# because of the time scale: the ram-pressure balance
# ρ_j(v_j − v_h)² = ρ_a v_h² gives a jet-head speed
# v_h = v_j/(1 + √(ρ_a/ρ_j)) = 800/(1+√0.1) ≈ 608, so at the paper's final time
# t = 2e-3 the head is near y ≈ 1.2 and has not reached the top of the 1.5-tall
# box — nothing has had time to reflect off any boundary and come back.
#
# INFLOW. The beam is injected at Mach 800, so EVERY characteristic of the MHD
# system enters the domain through the nozzle and the whole nine-component
# conserved state is prescribed:
#
#   ρ = γ,  v = (0, 800, 0),  p = 1,  B = (0, B_a, 0),  ψ = 0
#
# (aj_jet_state() in user_flux.jl). B is prescribed as well as the fluid state:
# the field is uniform and y-aligned in both the ambient medium and the beam,
# so injecting B = (0, B_a, 0) keeps ∇·B = ∂_y B_a = 0 exactly at the inlet,
# whereas leaving B free there would let the inlet become a source of
# divergence error that the GLM ψ then has to clean.
#
# THE NOZZLE LIP. |x| = 0.05 falls exactly on an element boundary in every mesh
# that ships with this case (h = 0.025 and h = 0.01 both divide 0.05), so there
# IS a node at |x| = 0.05 sitting on the discontinuity of the boundary datum.
# The closed test |x| <= 0.05 is used, i.e. that node is injected: the inflow
# patch is then exactly the closed segment of width 0.1 that the paper's
# nozzle is. The tolerance is absolute and 10⁻¹⁰, nine orders below the
# smallest element, because the lip coordinate comes out of the mesh as
# 0.050000000000000044 rather than 0.05.
#
# THE NOZZLE LIP IS WHERE THE FIRST RUN DIED, and the measurement is exact: the
# GLOBAL first realizability repair of the whole run was at
#
#     (x, y) = (-0.075, 0.0)   on RHS call 3
#
# i.e. on the bottom boundary, one element OUTSIDE the lip, during the first time
# step. Nothing propagated there — the fastest signal travels c_h*dt = 4e-4 in a
# step, 1.6 % of an element — so the boundary condition put it there, exactly as a
# 446 m/s shear held open by a Dirichlet condition put rampCaoEtAl2021_M7's first
# repair one LGL interval below its top boundary.
#
# The mechanism is the CLAMPED/FREE INTERFACE, not the lip coordinate itself. A
# top-hat datum imposed strongly means the nodes at |x| <= 0.05 are pinned to
# rhoE = 4.48e5 while the node next door is free at the ambient 102.5 — a jump of
# 4400x inside one spectral element, with the element's collocation derivative
# reading d(rhoE)/dx ~ 1.8e7 across it. The Gibbs response is of the jump's own
# size, so a node whose rhoE is 102.5 receives an excursion of O(1e5). It cannot
# survive that, and no amount of viscosity fixes it: DynSGS was measured at
# nu = 2.1013 = C_max*Delta*lambda, i.e. SATURATED AT ITS CAP, in the same run.
#
# THE FIX (JEXPRESSO_AJ_SMOOTH, on by default): blend the injected state into the
# ambient one with a COMPACTLY SUPPORTED profile centred on the lip,
#
#     phi(x) = 1                        |x| <= x0 - s
#            = 1 - smootherstep(xi)     x0 - s < |x| < x0 + s,  xi = (|x|-(x0-s))/2s
#            = 0                        |x| >= x0 + s
#
# with smootherstep(xi) = 6xi^5 - 15xi^4 + 10xi^3, and the Dirichlet patch widened
# to |x| <= x0 + s to match. Three properties, all of which the previous tanh
# version lacked:
#
#   * phi REACHES EXACTLY ZERO at the patch edge, so the outermost CLAMPED node
#     holds precisely the ambient state — the same thing its free neighbour holds.
#     There is no clamped/free jump anywhere. A tanh cannot do this: it only
#     decays, so phi at the patch edge is small but nonzero, and "small" is not
#     small enough here — 0.02 % of the beam's rhoE is 100, which is the ENTIRE
#     ambient rhoE. That is why tapering inside the nozzle only (the first
#     implementation) merely halved the jump and did not remove it.
#   * phi is C2 at both ends (smootherstep has zero first AND second derivative
#     there), so a spectral element can represent the datum without ringing at
#     the joins. A C0 kink inside an element is what cost rampCaoEtAl2021_M7 a
#     17.5x pressure amplification.
#   * it is ANTISYMMETRIC about |x| = x0, so the injected mass and momentum flux
#     are unchanged to second order: the beam keeps its width. (The alternative,
#     tapering inward to zero AT the lip, keeps the Dirichlet patch exactly equal
#     to the paper's nozzle but narrows the full-strength beam to |x| <= x0 - 2s.
#     Preserving the jet was judged the more important of the two.)
#
# THE COST, stated plainly: the Dirichlet patch is s wider than the paper's nozzle
# on each side, so a strip of the bottom boundary of that width is held at a blend
# of beam and INITIAL ambient instead of being open. That is a real modelling
# error where the cocoon later flows back down — it is 1.25 % of the bottom
# boundary per side at the default s — and it is much smaller than the error of
# imposing a datum the discretization cannot represent.
#
# s DEFAULTS TO ONE ELEMENT, resolved from the mesh by initialize.jl, so the
# transition 2s spans two elements and each element's P4 polynomial sees half of
# the smootherstep rather than all of it.
#
# NOTE WHAT s DOES NOT BUY. The worst node-to-node step in rhoE improves only from
# 2.1x to 2.3x better than the top hat's between s = h/2 and s = h, because
# rhoE ~ rho(phi)*(phi*u)^2 is intrinsically steep near phi = 1 whatever the width.
# The point of s is CONTINUITY and per-element resolvability, not a smaller
# gradient: the datum really does span 4400x across the nozzle shoulder in any
# representation, and a steep RESOLVED profile is what DynSGS is for, whereas a
# discontinuity between a clamped node and a free one is not. Beyond s ~ 2h the
# full-strength core |x| <= x0 - s vanishes and the beam stops being the paper's.
#
# And the injected flux is unchanged: the integral of phi over the boundary is x0
# to five digits for every s, because the profile is antisymmetric about the lip.
# So s alters the beam's shoulder shape, not how much mass or momentum enters.
# JEXPRESSO_AJ_SMOOTH=0 restores the paper's exact top hat — which is the faithful
# condition and is measured NOT to run; see README.md 10-11.
#---------------------------------------------------------------------------------

if !@isdefined(AJ_LIP_TOL)
    const AJ_LIP_TOL = 1.0e-10
end

# Half-width of the Dirichlet patch: the nozzle plus the smoothing reach, so that
# phi has already fallen to exactly 0 by the outermost clamped node.
@inline aj_patch_halfwidth() = AJ_XNOZZLE + max(aj_smooth[], 0.0)

# The blend factor. s <= 0 is the paper's top hat: 1 on the closed nozzle, 0 off it.
@inline function aj_nozzle_phi(x)
    s  = aj_smooth[]
    ax = abs(x)
    s <= 0.0 && return ax <= AJ_XNOZZLE + AJ_LIP_TOL ? 1.0 : 0.0
    ax <= AJ_XNOZZLE - s && return 1.0
    ax >= AJ_XNOZZLE + s && return 0.0
    ξ = (ax - (AJ_XNOZZLE - s))/(2.0*s)              # 0 .. 1 across the transition
    # 1 - smootherstep(ξ): C2 at both ends, exactly 1 and exactly 0 there
    return 1.0 - ξ*ξ*ξ*(10.0 + ξ*(-15.0 + 6.0*ξ))
end

# The conserved state injected at a nozzle node of abscissa x.
#
# The blend is applied to the PRIMITIVES (ρ, v, p, B) and ρE is rebuilt from them.
# Blending the conserved variables instead puts a spurious pressure spike on the
# lip: at φ = ½ the kinetic energy of the mean momentum is not the mean of the
# kinetic energies, and p would come out ≈ 8000 instead of 1.
@inline function aj_nozzle_state(x)
    Ba = aj_Ba[]
    φ  = aj_nozzle_phi(x)

    ρ = φ*AJ_RHO_JET + (1.0 - φ)*AJ_RHO_AMB
    v = φ*aj_ujet[]                                  # the ambient is at rest
    p = φ*AJ_P_JET   + (1.0 - φ)*AJ_P_AMB            # = 1 either way
    # Bx = Bz = ψ = 0 and By = B_a in both states, so the blend leaves B alone.
    ρE = p/(γ_mhd - 1.0) + 0.5*ρ*v*v + 0.5*Ba*Ba

    return (ρ, 0.0, ρ*v, ρE, 0.0, 0.0, Ba, 0.0, 0.0)
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "bottom" && abs(coords[1]) <= aj_patch_halfwidth() + AJ_LIP_TOL
        s = aj_nozzle_state(coords[1])
        for ieq = 1:9
            qbdy[ieq] = s[ieq]
        end
    end
    # Everything else — "left", "right", "top", and the bottom boundary outside
    # the nozzle — is outflow: impose nothing, leave qbdy at its sentinel.

    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,
                            qe, ::PERT)
    error(" problems/MHD/astroJetWuShu2018: PERT() solution variables are not supported.")
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, inputs)
    flux = zeros(size(q,2),1)
    return flux
end
