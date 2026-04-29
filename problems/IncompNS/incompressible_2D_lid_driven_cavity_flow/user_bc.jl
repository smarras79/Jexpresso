

# function user_bc_dirichlet!(Coef, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, iel, params, ip_bd)

#     #= xmin = params.mesh.xmin; 
#     xmax = params.mesh.xmax; 
#     ymax = params.mesh.ymax;
#     ymin = params.mesh.ymin;
#     X = params.mesh.x
#     Y = params.mesh.y
#     w = params.uaux
#     thom = 1
#     woods = 2

#     type = woods

#     if ((x == xmin || x == xmax) && (y != ymin && y != ymax))
#         ip = find_lr(params, x, y, iel)
#         temp = Coef[ip]
#         bd = Coef[ip_bd]

#         if(x == xmin)
#             dx = abs(xmin - X[ip])
#         else
#             dx = abs(xmax - X[ip])
#         end

#         if (type == thom)
#             qbdy[1] = 2*(temp - bd) / (dx^2)
#         else
#             qbdy[1] = -0.5*w[ip] - 3*(bd- temp)/(dx^2)
#         end

#     elseif ((x != xmin && x != xmax) && (y == ymin || y == ymax))
#         ip = find_ud(params, x, y, iel)
#         temp = Coef[ip]
#         bd = Coef[ip_bd]

#         if(y == ymin)
#             dy = abs(ymin - Y[ip])
            
#             if (type == thom)
#                 qbdy[1] = 2*(temp - bd) / (dy^2)
#             else
#                 qbdy[1] = -0.5*w[ip] - 3*(bd- temp)/(dy^2)
#             end

#         else
#             dy = abs(ymax - Y[ip])
#             v = velocity(x,xmin,xmax)

#             if (type == thom)
#                 qbdy[1] = 2*(temp - bd) / (dy^2) - 2*v / dy
#             else
#                 qbdy[1] = -0.5*w[ip] - 3*(bd- temp)/(dy^2) - 3*v/dy
#             end

#         end
        
#     else # corner
#         ip = find_ud(params, x, y, iel)
#         temp = Coef[ip]
#         bd = Coef[ip_bd]
       
#         if(y == ymin) # down
#             dy = abs(ymin - Y[ip])

#             if (type == thom)
#                 qbdy[1] = 2*(temp - bd) / (dy^2)
#             else
#                 qbdy[1] = -0.5*w[ip] - 3*(bd- temp)/(dy^2)
#             end

#         else # up
#             dy = abs(ymax - Y[ip])
#             v = velocity(x,xmin,xmax)

#             if (type == thom)
#                 qbdy[1] = 2*(temp - bd) / (dy^2) - 2*v / dy
#             else
#                 qbdy[1] = -0.5*w[ip] - 3*(bd- temp)/(dy^2) - 3*v/dy
#             end

#         end

#     end

#     if (t == 0.0)
#         qbdy[1] = 0.0
#     end =#
#     nothing
# end

# function user_bc_dirichlet!(q, x, y, t, bdy_edge_type, qbdy, nx_l, ny_l, qe,SOL_VARS_TYPE)
#     nothing

# end

# function find_lr(params, x_bd, y_bd, iel)

#     connijk = params.mesh.connijk
#     ngl = params.mesh.ngl
#     x = params.mesh.x
#     y = params.mesh.y
#     pos = zeros(2,ngl-1)
#     id = 1

#     for j = 1:ngl
#         for i = 1:ngl
#             ip = connijk[iel,i,j]
#             xx = x[ip]
#             yy = y[ip]
#             if(round(yy,digits = 5) == round(y_bd,digits = 5) && round(xx,digits = 5) != round(x_bd,digits = 5))
#                 pos[1,id] = ip
#                 pos[2,id] = abs(xx-x_bd)
#                 id = id+1
#             end
#         end
#     end

#     min_val, min_index = findmin(pos[2,:])

#     return round(Int64,pos[1,min_index])

# end

# function find_ud(params, x_bd, y_bd, iel)

#     connijk = params.mesh.connijk
#     ngl = params.mesh.ngl
#     x = params.mesh.x
#     y = params.mesh.y

#     pos = zeros(2,ngl-1)
#     id = 1

#     for j = 1:ngl
#         for i = 1:ngl
#             ip = connijk[iel,i,j]
#             xx = x[ip]
#             yy = y[ip]
#             if(round(yy,digits = 5) != round(y_bd,digits = 5) && round(xx,digits = 5) == round(x_bd,digits = 5))
#                 pos[1,id] = ip
#                 pos[2,id] = abs(yy-y_bd)
#                 id = id+1
#             end
#         end
#     end

#     min_val, min_index = findmin(pos[2,:])

#     return round(Int64,pos[1,min_index])

# end

# function velocity(x,xmin,xmax)
    
#     return 1
#     # return 16*(x-xmin)^2*(x-xmax)^2

# end
function user_bc_dirichlet!(x, y, t, tag, velocity)

    if (tag == "LID")
        velocity[1] = 1
        velocity[2] = 0
    else
        velocity[:] .= 0
    end

end