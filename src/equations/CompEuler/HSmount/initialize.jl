using Base

function initialize(SD::NSD_2D, PT::CompEuler, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

                """
    @info " Initialize fields for 2D CompEuler with θ equation ........................ "
    
    PhysConst = PhysicalConst{Float64}()
    
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=4)
    
    θref = 250.0 #K
    T0   = θref
    p0   = 100000.0
    
    N    = PhysConst.g/sqrt(PhysConst.cp*T0)
    N2   = N*N
    
    for iel_g = 1:mesh.nelem
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            ip = mesh.connijk[i,j,iel_g]
            y = mesh.y[ip]
            
            θ    = θref*exp(N*y/PhysConst.g)            
            p    = p0*(1.0 + PhysConst.g2*(exp(-y*N2/PhysConst.g) - 1.0)/(PhysConst.cp*θref*N2))^PhysConst.cpoverR
            ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p) #kg/m³
            ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

            u = 20.0
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

    for iel_g = 1:mesh.nelem_semi_inf
        for j=1:mesh.ngr, i=1:mesh.ngl

            ip = mesh.connijk_lag[i,j,iel_g]
            x, y = mesh.x[ip], mesh.y[ip]
            
            p    = p0*(1.0 + PhysConst.g2*(exp(-y*N2/PhysConst.g) - 1.0)/(PhysConst.cp*θref*N2))^PhysConst.cpoverR
            ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p) #kg/m³
            ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³
            
            u = 20.0
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
    
    @info "Initialize fields for system of 2D CompEuler with θ equation ........................ DONE"
    
    return q
end
