#---------------------------------------------------------------------------------
# Perforated plastic plate — Artioli, Beirão da Veiga, Lovadina & Sacco,
# Comput. Mech. 60 (2017), Sec. 4.3 (arXiv:1701.06676), after the classic
# perforated strip of Zienkiewicz & Taylor.
#
# A strip 2L x 2H = 200 x 360 mm with a central hole of diameter 2R = 100 mm,
# plane strain, J2 plasticity (paper Sec. 3.2) with
#
#     E = 7000 kg/mm²,   ν = 0.3,   σy0 = 24.3 kg/mm²,   no hardening.
#
# A quarter is modelled. The normal displacement is restrained on the two
# symmetry edges, on the lateral edge x = L and on the top edge y = H, where it
# is prescribed: the top is pulled up by δ = 2 mm in 400 equal increments
# (paper Fig. 5.3a). Reported: the horizontal displacement of A = (R, 0) and
# the vertical displacement of B = (0, R) at the end (Table 1), the Newton
# iterations per increment (Table 2) and the reaction sum on the top edge
# (Fig. 5.4; the paper's force is that of a 10 mm thick plate, 10x the per-mm
# force this plane-strain model computes).
#
# Runs on the quasi-static solid path with the history at the quadrature
# points (src/kernel/solvers/quasistatic_solid.jl, model = :j2_plasticity).
# t ∈ [0, 1] is a load parameter: J2 plasticity is rate independent.
#
# A METHOD and not `const`s, as in the sibling cases (src/run.jl includes a
# case into the Jexpresso module).
#---------------------------------------------------------------------------------
@inline function perforated_plate_properties()

    E    = 7000.0     # Young's modulus      [kg/mm²]
    ν    = 0.3        # Poisson's ratio
    σy0  = 24.3       # yield stress          [kg/mm²]
    Hiso = 0.0        # no hardening: perfectly plastic
    Hkin = 0.0

    L    = 100.0      # half width   [mm]   (must match the mesh)
    H    = 180.0      # half length  [mm]
    R    = 50.0       # hole radius  [mm]
    δmax = 2.0        # displacement of the top edge at t = 1 [mm]

    K = E/(3.0*(1.0 - 2.0*ν))
    G = E/(2.0*(1.0 + ν))

    return (; E, ν, σy0, Hiso, Hkin, L, H, R, δmax, K, G)
end

#
# Hook read by the quasi-static solver.
#
function user_solid_material()
    P = perforated_plate_properties()
    return (; model = :j2_plasticity, E = P.E, ν = P.ν, σy0 = P.σy0, Hiso = P.Hiso, Hkin = P.Hkin)
end

# The displacement prescribed on the top edge, ramped linearly with t ∈ [0, 1].
@inline top_displacement(t) = perforated_plate_properties().δmax*t
