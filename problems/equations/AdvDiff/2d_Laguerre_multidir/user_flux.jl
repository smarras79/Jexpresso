function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL,ip; neqs=1)
    
    u = 2.0
    v = 0.0
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
    
end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::PERT,ip; neqs=1)
   
    x=mesh.x[ip]
    y=mesh.y[ip]
    r =sqrt((x-5.0)^2+(y-5.0)^2)
    u = 0.0
    v = -0.0
    if (r > 1.0 && x >6.0)
      u = 2.0
      v = 0.0
    elseif (r > 1.0 && x <4.0)
      u = -2.0
      v = 0.0
    elseif (r > 1.0 && y >6.0)
      u = 0.0
      v = 2.0
    elseif (r > 1.0 && y <4.0)
      u = 0.0
      v = -2.0
    elseif (r <= 1.0)
      u = 2*(x-5.0)
      v = 2*(y-5.0)
    end   
   
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
    
end


function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::NCL, ::AbstractPert; neqs=1)
    
    u = 0.5
    v = 1.0
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
end
