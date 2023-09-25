using Base

function initialize(SD::NSD_2D, PT::CompEuler, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    
    @info " Initialize fields for 2D CompEuler with θ equation ........................ "

    PhysConst = PhysicalConst{Float64}()
    
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=4)

     if (inputs[:case] === "mountain")
       
         @info "   main domain height, and laguerre height", mesh.ymax, maximum(mesh.y)
  
         θ0 = 250.0 #K
         T0 = θ0
         p0 = 100000.0
         
         N    = 0.01
         N2   = N*N
         g2 = PhysConst.g*PhysConst.g
         
         for iel_g = 1:mesh.nelem
             for j=1:mesh.ngl, i=1:mesh.ngl
                 
                 ip = mesh.connijk[iel_g,i,j]
                 y  = mesh.y[ip]
                 
                 auxi = PhysConst.Rair*θ0
                 p    = p0*exp(-PhysConst.g*y/auxi)
                 θ    = θ0*exp(N2*y/PhysConst.g)
                 
                 ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p) #kg/m³
                 ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p) #kg/m³
                 
                 u = 0.0
                 v = 0.0
                 
                 q.qn[ip,1]   = ρ
                 q.qn[ip,2]   = ρ*u
                 q.qn[ip,3]   = ρ*v
                 q.qn[ip,4]   = ρ*θ
                 q.qn[ip,end] = p   #pressure
                 
                 #@info θ, ρ, p, y, mesh.x[ip]
                 #Store initial background state for plotting and analysis of pertuebations
                 q.qe[ip,1] = ρref
                 q.qe[ip,2] = ρref*u
                 q.qe[ip,3] = ρref*v
                 q.qe[ip,4] = ρref*0          
                 q.qe[ip,end] = p   #pressure
                 
             end
         end
         
     elseif(inputs[:case] === "rtb")
         
        xc = (maximum(mesh.x) + minimum(mesh.x))/2
        yc = 2500.0 #m
        r0 = 2000.0 #m
        
        θref = 300.0 #K
        θc   =  0.0# 2.0 #K
        for iel_g = 1:mesh.nelem
            for j=1:mesh.ngl, i=1:mesh.ngl
                
                ip = mesh.connijk[iel_g,i,j]
                x, y = mesh.x[ip], mesh.y[ip]
                r = sqrt( (x - xc)^2 + (y - yc)^2 )
                
                Δθ = 0.0 #K
                if r < r0
                    Δθ = θc*(1.0 - r/r0)
                end
                θ = θref + Δθ
                p    = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
                pref = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p) #kg/m³
                ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

                u = 0.0
                v = 0.0
                
                q.qn[ip,1] = ρ
                q.qn[ip,2] = ρ*u
                q.qn[ip,3] = ρ*v
                q.qn[ip,4] = ρ*θ
                q.qn[ip,end] = p

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1] = ρref
                q.qe[ip,2] = u
                q.qe[ip,3] = v
                q.qe[ip,4] = ρref*θref
                q.qe[ip,end] = pref   #pressure
            end
        end
        
     end
    @info "Initialize fields for system of 2D CompEuler with θ equation ........................ DONE"
    
    return q
end
