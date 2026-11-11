#---------------------------------------------------------------------------------
# Body force per unit volume, read by the quasi-static solver at every
# quadrature point: user_solid_body_force(x, y, t) -> (bx, by).
#
# None in this problem: the cylinder is loaded only by the pressure on its
# inner surface (user_bc.jl). Self-weight would be (0.0, -ρ*g).
#---------------------------------------------------------------------------------
user_solid_body_force(x, y, t) = (0.0, 0.0)
