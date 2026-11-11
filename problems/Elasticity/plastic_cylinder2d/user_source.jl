#---------------------------------------------------------------------------------
# Body force per unit volume: none. The cylinder is loaded only through the
# displacement of its bore (user_bc.jl).
#---------------------------------------------------------------------------------
user_solid_body_force(x, y, t) = (0.0, 0.0)
