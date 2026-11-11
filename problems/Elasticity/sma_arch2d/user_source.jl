#---------------------------------------------------------------------------------
# Body force per unit volume: none. The arch is loaded only through the
# traction on its free end (user_bc.jl).
#---------------------------------------------------------------------------------
user_solid_body_force(x, y, t) = (0.0, 0.0)
