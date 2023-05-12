#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("../../io/mod_inputs.jl")
include("../mesh/metric_terms.jl")
include("../mesh/mesh.jl")
include("../operators/operators.jl")
include("../abstractTypes.jl")
include("../bases/basis_structs.jl")
include("../infrastructure/element_matrices.jl")
include("../infrastructure/Kopriva_functions.jl")
include("../infrastructure/2D_3D_structures.jl")
include("custom_bcs.jl")

function apply_periodicity!(SD::NSD_1D, rhs, qp, mesh, inputs, QT, metrics, ψ, dψ, ω, t, nvars)

    #NOTICE: " apply_periodicity!() in 1D is now only working for nvars=1!"
    
    #
    # 1D periodic
    #
    qp[mesh.npoin_linear,:] .= 0.5*(qp[mesh.npoin_linear,:] .+ qp[1,:])
    qp[1,:] = qp[mesh.npoin_linear,:]

end

function apply_boundary_conditions!(SD::NSD_1D, rhs, qp, mesh,inputs, QT, metrics, ψ, dψ, ω, t, nvars;L=zeros(1,1))
    #If Neumann conditions are needed compute gradient
    calc_grad = false
    #   for key in keys(inputs)
    #     if (inputs[key] == "dirichlet" || inputs[key] == "neumann" || inputs[key] == "dirichlet/neumann")
    calc_grad = true
    #    end
    #  end
    if (inputs[:lneumann])        
        gradq = zeros(Float64, mesh.npoin,nvars)
    else
        gradq = zeros(Float64, 1,1)
    end
    #TODO remake build custom_bcs for new boundary data
    #if (calc_grad)
    #    gradq = build_gradient(SD, QT::Inexact, qp, ψ, dψ, ω, mesh, metrics,gradq,nvars)
    build_custom_bcs!(t,mesh,qp,gradq,rhs,SD,nvars,metrics,ω,dirichlet!,neumann,L,inputs)
    #end
end

function apply_boundary_conditions!(SD::NSD_2D, rhs, qp, mesh, inputs, QT, metrics, ψ, dψ, ω, t, nvars;L=zeros(1,1))
    #If Neumann conditions are needed compute gradient
    #calc_grad = false
    #   for key in keys(inputs)
    #     if (inputs[key] == "dirichlet" || inputs[key] == "neumann" || inputs[key] == "dirichlet/neumann")
    #calc_grad = true
    #    end
    #  end
    #nface = size(mesh.bdy_edge_comp,1)
    #dqdx_st = zeros(nvars,2)
    #q_st = zeros(nvars,1)
    if (inputs[:lneumann])
        gradq = zeros(Float64, 2, mesh.npoin, nvars)
    else
        gradq = zeros(Float64, 1, 1, 1)
    end
    #flux_q = zeros(mesh.ngl,nface,2,nvars)
    #exact = zeros(mesh.ngl,nface,nvars)
    #penalty =0.0#50000
    #nx = metrics.nx
    #ny = metrics.ny
    ##TODO remake build custom_bcs for new boundary data
    ##if (calc_grad)
    ##    gradq = build_gradient(SD, QT::Inexact, qp, ψ, dψ, ω, mesh, metrics,gradq,nvars)
    build_custom_bcs!(t,mesh,qp,gradq,rhs,SD,nvars,metrics,ω,dirichlet!,neumann,L,inputs)
    #end
    
end

function build_custom_bcs!(t,mesh,q,gradq,rhs,::NSD_1D,nvars,metrics,ω,dirichlet!,neumann,L,inputs)

    for ip in [1, mesh.npoin_linear]
        x = mesh.x[ip]
        if (ip == 1)
            k=1
            iel = 1
        else
            k=mesh.ngl
            iel = mesh.nelem
        end
        qbdy = zeros(size(q,2),1)
        qbdy[:] .= 4325789.0        
        flux = zeros(Float64, size(q,2),1)
        if (inputs[:luser_bc])
            qbdy = dirichlet!(q[ip,:], x,t,mesh,metrics,qbdy)
            if (inputs[:lneumann])
                flux = (ω[k]*neumann(q[ip,:],gradq[ip,:],x,t,mesh,metrics))
            end
        else
            q[ip,:] .= 0.0
        end
        
        rhs[k,iel,:] .= rhs[k,iel,:] .+ flux[:]
        for var =1:size(q,2)
            if !(AlmostEqual(qbdy[var],4325789.0))
                rhs[k,iel,var] = 0.0
                q[ip,var] = qbdy[var]
            end
        end
        if (size(L,1)>1)
            for ii=1:mesh.npoin
                L[ip,ii] = 0.0
            end
            L[ip,ip] =1.0
        end
    end
end

function build_custom_bcs!(t,mesh,q,gradq,rhs,::NSD_2D,nvars,metrics,ω,dirichlet!,neumann,L,inputs)

    
    for iedge = 1:size(mesh.bdy_edge_in_elem,1)
        iel = mesh.bdy_edge_in_elem[iedge]
        comp = mesh.bdy_edge_comp[iedge]
        for k=1:mesh.ngl
            if (mesh.bdy_edge_type[iedge] != "periodic1" && mesh.bdy_edge_type[iedge] !="periodic2")
                tag = mesh.bdy_edge_type[iedge]
                ip = mesh.poin_in_bdy_edge[iedge,k]
                mm=1
                ll=1
                for ii=1:mesh.ngl
                    for jj=1:mesh.ngl
                        if (mesh.connijk[ii,jj,iel] == ip)
                            m=jj
                            l=ii
                        end
                    end
                end
                flux = zeros(size(q,2),1)
                x = mesh.x[ip]
                y = mesh.y[ip]
                qbdy = zeros(size(q,2),1)
                qbdy[:] .= 4325789.0 
                #flags = zeros(size(q,2),1)
                flux = zeros(Float64, size(q,2),1)
                if (inputs[:luser_bc])
                    #q[ip,:], flags = dirichlet!(q[ip,:],x,y,t,mesh,metrics,tag,qbdy)
                    qbdy = dirichlet!(q[ip,:],x,y,t,mesh,metrics,tag,qbdy)
                    if (inputs[:lneumann])
                        flux .= (ω[k]*metrics.Jef[k,iedge]).*neumann(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag)
                    end
                else
                    q[ip,:] .= 0.0
                end
                rhs[ll,mm,iel,:] .= rhs[ll,mm,iel,:] .+ flux[:]
                
                for var =1:size(q,2)
                    if !(AlmostEqual(qbdy[var],4325789.0))
                        #@info var,x,y,qbdy[var]
                        rhs[ll,mm,iel,var] = 0.0
                        q[ip,var] = qbdy[var]
                    end
                end
                if (size(L,1)>1)
                    for ii=1:mesh.npoin
                        L[ip,ii] = 0.0
                    end
                    L[ip,ip] =1.0
                end
            end
        end
    end
end


