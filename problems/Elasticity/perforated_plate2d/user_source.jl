#---------------------------------------------------------------------------------
# Body force per unit volume: none. The plate is loaded only through the
# displacement of its top edge (user_bc.jl).
#---------------------------------------------------------------------------------
user_solid_body_force(x, y, t) = (0.0, 0.0)
