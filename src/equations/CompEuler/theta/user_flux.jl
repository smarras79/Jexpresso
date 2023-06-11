function user_flux(SD::NSD_2D, q::Array, mesh::St_mesh; neqs=4)
    PhysConst = PhysicalConst{Float64}()
    
    F = zeros(neqs)
    G = zeros(neqs)
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    
    F[1] = ρu
    F[2] = ρu*ρu/ρ + Press
    F[3] = ρu*ρv/ρ
    F[4] = ρu*ρθ/ρ

    F[1] = ρv
    F[2] = ρv*ρu/ρ
    F[3] = ρv*ρv/ρ + Press
    F[4] = ρv*ρθ/ρ
    
    return F, G
end


function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::SubArray{Float64}, mesh::St_mesh; neqs=4)
    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ

    Press = 0.0   
    perfectGasLaw_ρθtoP!(Press, PhysConst, ρ=ρ, θ=θ)
    
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρv*u
    F[4] = ρθ*u

    F[1] = ρv
    F[2] = ρu*v
    F[3] = ρv*v + Press
    F[4] = ρθ*v
    
end



primitives(x,y) = x/y
primitives2flux(x,y,z) = x.*y .+ z
function user_flux!(F::SubArray{Float64}, G::SubArray{Float64},
                    SD::NSD_2D,
                    ρ, ρu, ρv, ρθ,
                    mesh::St_mesh; neqs=4)
#function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::Array{Float64}, mesh::St_mesh; neqs=4)
    PhysConst = PhysicalConst{Float64}()
    
    Press = 0.0
    perfectGasLaw_ρθtoP!(Press, PhysConst, ρ=ρ, θ=ρθ/ρ)
  #=  a=q[2].+Press
    Fflux = @MVector [q[2]
                      a
                      q[2]
                      q[2]]
    
    b=q[3].+Press
    Gflux =  @MVector [q[3]
                       q[3]
                       b
                       q[3]]
   =# 
 #=   Fflux = @MVector [q[2]
                      q[2].*q[2]./q[1] .+ Press
                      q[3].*q[2]./q[1]
                      q[4].*q[2]./q[1]]
    
    Gflux =  @MVector [q[3]
                       q[3].*q[2]./q[1]
                       q[3].*q[3]./q[1] .+ Press
                       q[4].*q[3]./q[1]]
    =#
    #Fflux = @MVector zeros(4)
    #Gflux = @MVector zeros(4)
    #prim = primitive(SA[ρu, ρv, ρθ], ρ)

    #Fflux = primitives2flux(SA[ρu, ρv, ρθ], SA[ρu, ρv, ρθ], SA[Press, 0.0, 0.0])
    #Gflux = primitives2flux(SA[ρu, ρv, ρθ], SA[ρu, ρv, ρθ], SA[0.0, Press, 0.0])

    F = primitives2flux(SA[1.0, 1.0, 1.0, 1.0], SA[1.0, 1.0, 1.0, 1.0], SA[0.0, 1.0, 0.0, 0.0])
    G = primitives2flux(SA[1.0, 1.0, 1.0, 1.0], SA[1.0, 1.0, 1.0, 1.0], SA[0.0, 0.0, 1.0, 0.0])
        
    #primitive2flux!(F, G, q)

#    F[1] = ρu, F[2:end] .= Fflux[1:end]
#    G[1] = ρv, G[2:end] .= Gflux[1:end]
    
    #F[2] .= ρu*u + Press
    #F[3] .= ρv*u
    #F[4] .= ρθ*u
    
    #G[1] .= ρv
    #G[2] .= ρu*v
    #G[3] .= ρv*v + Press
    #G[4] .= ρθ*v
    
    #return Fflux, Gflux
end
function newuser_flux!(F::Array{Float64}, G::Array{Float64}, SD::NSD_2D, q::Array, iel, mesh::St_mesh; neqs=4)
    
    PhysConst = PhysicalConst{Float64}()

    for j=1:mesh.ngl, i=1:mesh.ngl
        
        ip = mesh.connijk[i,j,iel]
        
        ρ  = q[ip,1]
        ρu = q[ip,2]
        ρv = q[ip,3]
        ρθ = q[ip,4]
        θ  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        
        Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
        
        F[i,j,1] = ρu
        F[i,j,2] = ρu*u + Press
        F[i,j,3] = ρv*u
        F[i,j,4] = ρθ*u
        
        G[i,j,1] = ρv
        G[i,j,2] = ρu*v
        G[i,j,3] = ρv*v + Press
        G[i,j,4] = ρθ*v
        
    end
end
