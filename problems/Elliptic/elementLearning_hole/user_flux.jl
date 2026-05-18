function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4)

    PhysConst = PhysicalConst{Float64}()

    F[1] = 0.0
    F[2] = 0.0
    F[3] = 0.0
    F[4] = 0.0

    G[1] = 0.0 
    G[2] = 0.0
    G[3] = 0.0
    G[4] = 0.0    
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4)

    PhysConst = PhysicalConst{Float64}()
    
    F[1] = 0.0
    F[2] = 0.0
    F[3] = 0.0
    F[4] = 0.0

    G[1] = 0.0 
    G[2] = 0.0
    G[3] = 0.0
    G[4] = 0.0
end


function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::NCL, ::AbstractPert; neqs=4)
    
    PhysConst = PhysicalConst{Float64}()
                
    F[1] = 0.0
    F[2] = 0.0
    F[3] = 0.0
    F[4] = 0.0

    G[1] = 0.0 
    G[2] = 0.0
    G[3] = 0.0
    G[4] = 0.0

end
