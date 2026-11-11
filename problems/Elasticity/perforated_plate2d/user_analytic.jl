#---------------------------------------------------------------------------------
# References for the perforated plate. There is no closed form, so the solver
# compares nothing at run time (user_solid_exact returns `nothing`); the case
# is compared with the numbers Artioli et al. (2017) report, here.
#
#   Table 1   horizontal displacement of A = (R, 0) and vertical displacement
#             of B = (0, R) at the end of the loading [mm]. "Ref" is their
#             reference, "an overkilling discretization of u-p-εᵛ quadrilateral
#             mixed finite elements".
#   Table 2   average number of Newton iterations per load increment.
#---------------------------------------------------------------------------------

user_solid_exact(x, y, t) = nothing        # no closed form

function perforated_plate_paper_reference()
    return (; uA = 2.741, vB = 1.859,                                     # Table 1, Ref
              quadratic = (; VEM_Quad = (2.720, 1.851), VEM_Tri = (2.719, 1.852),
                             VEM_Voronoi = (2.714, 1.852), FEM_Q9 = (2.733, 1.853),
                             FEM_T6 = (2.719, 1.853)),
              linear    = (; VEM_Quad = (2.538, 1.819), VEM_Tri = (2.485, 1.823),
                             VEM_Voronoi = (2.540, 1.815), FEM_Q4 = (2.692, 1.836)),
              newton_quadratic = (; VEM_Quad = 6.43, VEM_Tri = 6.11, VEM_Voronoi = 6.21,
                                    FEM_Q9 = 6.70, FEM_T6 = 6.11))       # Table 2
end
