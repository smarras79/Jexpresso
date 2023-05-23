using Base

function initialize(SD::NSD_2D, PT::CompEuler, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

    """
    @info " Initialize fields for 2D CompEuler with θ equation ........................ "
    
    PhysConst = PhysicalConst{Float64}()
    
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=4)
    
    if (inputs[:case] === "rtb")

        xc = (maximum(mesh.x) + minimum(mesh.x))/2
        yc = 2500.0 #m
        r0   = 2000.0 #m
        
        θref = 300.0 #K
        θc   =   2.0 #K
        for iel_g = 1:mesh.nelem
            for j=1:mesh.ngl, i=1:mesh.ngl
                
                ip = mesh.connijk[i,j,iel_g]
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

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1] = ρref
                q.qe[ip,2] = u
                q.qe[ip,3] = v
                q.qe[ip,4] = ρref*θref
            end
        end
        
    else
        error(" ERROR: CompEuler: initialize.jl:\n assign value to inputs[:case]")
    end
    

    @info "Initialize fields for system of 2D CompEuler with θ equation ........................ DONE"

    return q
end
