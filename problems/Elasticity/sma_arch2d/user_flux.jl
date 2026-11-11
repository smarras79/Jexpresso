#---------------------------------------------------------------------------------
# Shape-memory alloy arch — Artioli, Beirão da Veiga, Lovadina & Sacco,
# Comput. Mech. 60 (2017), Sec. 4.4 (arXiv:1701.06676), Fig. 5.5 and Table 3;
# the device and the material of Evangelista, Marfia & Sacco (2009).
#
# A semicircular arch 3.5 ≤ r ≤ 4.5 mm, 0 ≤ θ ≤ π, clamped at θ = π (y = 0,
# x < 0), its other end (y = 0, 3.5 ≤ x ≤ 4.5) loaded by a uniform traction q
# along +x. Plane stress. The Souza–Auricchio–Petrini model (paper Sec. 3.3;
# src/kernel/solvers/quasistatic_solid_sma.jl) with
#
#     E = 53000 MPa, ν = 0.36, εL = 0.04, Mf = 223 K, h = 1000 MPa,
#     β = 2.1 MPa/K, σy0 = 50 MPa,
#
# under the load–temperature history of Table 3, five branches of the load
# parameter t:
#
#     t     0     1     2     3     4     5
#     q     0   qmax    0  −qmax    0     0         qmax = 60 N/mm
#     T   Troom Troom Troom Troom Troom Troom+80    Troom = 223 K = Mf
#
# piecewise linear in between, 40 equal increments in each of the first four
# branches and 10 in the heating one. At T = Mf a load leaves a transformation
# strain behind (h εL < R: no reverse transformation on unloading), the reversed
# load re-orients it, and heating by 80 K at zero load (τ = β·80 = 168 MPa > R)
# recovers the shape: the arch returns to its initial configuration.
#
# q is a load per unit length of the loaded end: in this unit-thickness 2D model
# a traction of q [N/mm²] on the 1 mm wide end face, the resultant q [N/mm]. The
# displacement reported (Fig. 5.6) is the horizontal one of A = (4.5, 0), the
# outer corner of the loaded end.
#
# Methods, not `const`s (src/run.jl includes a case into the Jexpresso module).
#---------------------------------------------------------------------------------
_sma_env(key, default) = String(strip(get(ENV, key, default)))
sma_verify_mode() = lowercase(_sma_env("JEXPRESSO_SMA_VERIFY", "")) == "elastic"

@inline function sma_arch_properties()

    E    = something(tryparse(Float64, _sma_env("JEXPRESSO_SMA_E", "")), 53000.0)   # Young's modulus [MPa]
    ν    = 0.36       # Poisson's ratio
    εL   = 0.04       # saturation transformation strain
    Mf   = 223.0      # martensite finish temperature      [K]
    h    = 1000.0     # transformation hardening           [MPa]
    β    = 2.1        # dτ/dT                              [MPa/K]
    σy0  = 50.0       # radius of the elastic domain is √(2/3) σy0 [MPa]

    a    = 3.5        # inner radius   [mm]   (must match the mesh)
    b    = 4.5        # outer radius   [mm]
    qmax = 60.0       # load amplitude [N/mm]
    # out-of-plane thickness the load q [N/mm] is spread over: the traction is
    # q/thick [N/mm²]. 1 mm: q is the resultant per unit thickness, as stated.
    thick = something(tryparse(Float64, _sma_env("JEXPRESSO_SMA_THICKNESS", "")), 1.0)
    Troom = 223.0     # [K]
    ΔT    = 80.0      # heating of the last branch [K]

    return (; E, ν, εL, Mf, h, β, σy0, a, b, qmax, thick, Troom, ΔT)
end

# Table 3: the load and the temperature, piecewise linear in t ∈ [0, 5].
function sma_arch_load(t)
    P = sma_arch_properties()
    t <= 1 && return  P.qmax*t
    t <= 2 && return  P.qmax*(2 - t)
    t <= 3 && return -P.qmax*(t - 2)
    t <= 4 && return -P.qmax*(4 - t)
    return 0.0
end

function sma_arch_temperature(t)
    P = sma_arch_properties()
    return t <= 4 ? P.Troom : P.Troom + P.ΔT*min(t - 4, 1.0)
end

# The hypothesis: plane stress (the paper), or plane strain for comparison.
sma_arch_hypothesis() = Symbol(lowercase(_sma_env("JEXPRESSO_SMA_HYPOTHESIS", "plane_stress")))

#
# Hook read by the quasi-static solver. In the verification mode
# (JEXPRESSO_SMA_VERIFY=elastic) the elastic domain is made so large that
# nothing transforms: the arch is linear elastic and has an exact solution
# (user_analytic.jl).
#
function user_solid_material()
    P = sma_arch_properties()
    return (; model = :sma, E = P.E, ν = P.ν, h = P.h, β = P.β, Mf = P.Mf, εL = P.εL,
              σy0 = sma_verify_mode() ? 1.0e12 : P.σy0,
              hypothesis = sma_arch_hypothesis(), temperature = sma_arch_temperature)
end
