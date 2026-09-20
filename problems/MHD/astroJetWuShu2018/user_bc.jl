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
# JEXPRESSO_AJ_SMOOTH=w > 0 blends the injected state into the ambient one
# across the lip with φ(x) = ½(1 − tanh((|x| − 0.05)/w)) — the paper's
# condition is the sharp top hat, w = 0, which is what runs by default. The
# blend is applied to the PRIMITIVE variables (ρ, v, p, B) and ρE is rebuilt
# from them: blending the conserved variables instead would put a spurious
# pressure spike on the lip, because at φ = ½ the kinetic energy of the mean
# momentum is not the mean of the kinetic energies (½(ρv)²/ρ evaluated on the
# blend is 8·10³ out of 2·10⁵ short of the blend of the two, i.e. p would come
# out ~8000 instead of 1).
#---------------------------------------------------------------------------------

if !@isdefined(AJ_LIP_TOL)
    const AJ_LIP_TOL = 1.0e-10
end

# The conserved state injected at a nozzle node of abscissa x.
@inline function aj_nozzle_state(x)
    Ba = aj_Ba[]
    w  = aj_smooth[]
    φ  = (w > 0.0) ? 0.5*(1.0 - tanh((abs(x) - AJ_XNOZZLE)/w)) : 1.0

    ρ = φ*AJ_RHO_JET + (1.0 - φ)*AJ_RHO_AMB
    v = φ*aj_ujet[]                                  # the ambient is at rest
    p = φ*AJ_P_JET   + (1.0 - φ)*AJ_P_AMB            # = 1 either way
    # Bx = Bz = ψ = 0 and By = B_a in both states, so the blend leaves B alone.
    ρE = p/(γ_mhd - 1.0) + 0.5*ρ*v*v + 0.5*Ba*Ba

    return (ρ, 0.0, ρ*v, ρE, 0.0, 0.0, Ba, 0.0, 0.0)
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    if tag == "bottom" && abs(coords[1]) <= AJ_XNOZZLE + AJ_LIP_TOL
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
