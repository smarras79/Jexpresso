using Base

function initialize(SD::NSD_2D, PT::CompEuler, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

                """
    @info " Initialize fields for 2D CompEuler with θ equation ........................ "
    

    PhysConst = PhysicalConst{Float64}()
    
    qvars = ("dρ", "dρu", "dρv", "dρθ")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    
    θref = 280.0 #K
    θ0 = 280.0
    T0   = θ0
    p0   = 100000.0
    
    N    = 0.01#PhysConst.g/sqrt(PhysConst.cp*T0)
    N2   = N*N
    
    for iel_g = 1:mesh.nelem
        for j=1:mesh.ngl, i=1:mesh.ngl
            
            ip = mesh.connijk[iel_g,i,j]
            y = mesh.y[ip]
            θ    = θref*exp(N2*y/PhysConst.g)         
            #if (y > 0.1) 
              p    = p0*(1.0 + PhysConst.g2*(exp(-y*N2/PhysConst.g) - 1.0)/(PhysConst.cp*θref*N2))^PhysConst.cpoverR
            #else
             # p = p0
            #end
            ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p) #kg/m³
            ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p) #kg/m³
            #=
            auxi = PhysConst.Rair*θ0
                 p    = p0*exp(-PhysConst.g*y/auxi)
                 θ    = θ0*exp(N2*y/PhysConst.g)
                 
                 ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p) #kg/m³
                 ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)
            =#
            u = 10.0
            v = 0.0
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[ip,1] = ρ-ρref
                q.qn[ip,2] = ρ*u - ρref*u
                q.qn[ip,3] = ρ*v - ρref*v
                q.qn[ip,4] = ρ*θ - ρref*θ
                q.qn[ip,end] = p
            else
                q.qn[ip,1] = ρ
                q.qn[ip,2] = ρ*u
                q.qn[ip,3] = ρ*v
                q.qn[ip,4] = ρ*θ
                q.qn[ip,end] = p
            end
            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,1] = ρref
            q.qe[ip,2] = ρref*u
            q.qe[ip,3] = ρref*v
            q.qe[ip,4] = ρref*θ
            q.qe[ip,end] = p
        end
    end

    #=for iel_g = 1:mesh.nelem_semi_inf
        for j=1:mesh.ngr, i=1:mesh.ngl

            ip = mesh.connijk_lag[iel_g,i,j]
            y = mesh.y[ip]
            θ    = θref*exp(N2*y/PhysConst.g)
            #if (y > 0.1)
              p    = p0*(1.0 + PhysConst.g2*(exp(-y*N2/PhysConst.g) - 1.0)/(PhysConst.cp*θref*N2))^PhysConst.cpoverR
            #else
             # p = p0
            #end
            ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p) #kg/m³
            ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p) #kg/m³
            #=
            auxi = PhysConst.Rair*θ0
                 p    = p0*exp(-PhysConst.g*y/auxi)
                 θ    = θ0*exp(N2*y/PhysConst.g)

                 ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p) #kg/m³
                 ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)
            =#
            u = 10.0
            v = 0.0
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[ip,1] = ρ-ρref
                q.qn[ip,2] = ρ*u - ρref*u
                q.qn[ip,3] = ρ*v - ρref*v
                q.qn[ip,4] = ρ*θ - ρref*θ
                q.qn[ip,end] = p
            else
                q.qn[ip,1] = ρ
                q.qn[ip,2] = ρ*u
                q.qn[ip,3] = ρ*v
                q.qn[ip,4] = ρ*θ
                q.qn[ip,end] = p
            end
            #Store initial background state for plotting and analysis of pertuebations
            q.qe[ip,1] = ρref
            q.qe[ip,2] = ρref*u
            q.qe[ip,3] = ρref*v
            q.qe[ip,4] = ρref*θ
            q.qe[ip,end] = p
        end
    end=#
    
    @info "Initialize fields for system of 2D CompEuler with θ equation ........................ DONE"
    
    return q
end
