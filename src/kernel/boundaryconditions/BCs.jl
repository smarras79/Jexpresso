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

function apply_periodicity!(SD::NSD_1D, qp, mesh, inputs)
    
    if (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="periodic" || haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="periodic")
        #
        # 1D periodic
        #
        qp[mesh.npoin_linear,:] .= 0.5*(qp[mesh.npoin_linear,:] .+ qp[1,:])
        qp[1,:] .= qp[mesh.npoin_linear,:]     
    elseif (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="dirichlet" || haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="dirichlet")
        #
        # Dirichlet q(1,t) = q(mesh.npoin_linear,t) = 0.0
        #
        qp[1,:] = 0.0
        qp[mesh.npoin_linear,:] = 0.0
    else
        #
        # 1D periodic
        #
        qp[mesh.npoin_linear,:] .= 0.5*(qp[mesh.npoin_linear,:] .+ qp[1,:])
        qp[1,:] .= qp[mesh.npoin_linear,:]
    end
end

function apply_boundary_conditions!(SD::NSD_1D, rhs, qp, mesh,inputs, QT, metrics, ψ, dψ, ω, t, nvars)
    nothing
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
    #flux_q = zeros(mesh.ngl,nface,2,nvars)
    #exact = zeros(mesh.ngl,nface,nvars)
    #penalty =0.0#50000
    #nx = metrics.nx
    #ny = metrics.ny
    #TODO remake build custom_bcs for new boundary data
    #if (calc_grad)
    gradq = zeros(2,mesh.npoin,nvars)
    #    gradq = build_gradient(SD, QT::Inexact, qp, ψ, dψ, ω, mesh, metrics,gradq,nvars)
    build_custom_bcs!(t,mesh,qp,gradq,rhs,SD,nvars,metrics,ω,dirichlet!,neumann,L,inputs)
    #end
end


function build_custom_bcs!(t,mesh,q,gradq,rhs,::NSD_2D,nvars,metrics,ω,dirichlet!,neumann,L,inputs)

    if (inputs[:luser_bc])
        for iedge = 1:size(mesh.bdy_edge_comp,1)
            iel = mesh.bdy_edge_in_elem[iedge]
            comp = mesh.bdy_edge_comp[iedge]
            if (mesh.bdy_edge_type[iedge] != "periodic1" && mesh.bdy_edge_type[iedge] !="periodic2")
                for k=1:mesh.ngl
                    tag = mesh.bdy_edge_type[iedge]
                    ip = mesh.poin_in_bdy_edge[iedge,k]
                    m=1
                    l=1
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            if (mesh.connijk[ii,jj,iel] == ip)
                                mm=jj
                                ll=ii
                            end
                        end
                    end
                    x = mesh.x[ip]
                    y = mesh.y[ip]
                    q[ip,:] = dirichlet!(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag)
                    flux = (ω[k]*metrics.Jef[k,iedge]).*neumann(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag)
                    rhs[l,m,iel,:] .= rhs[l,m,iel,:] .+ flux[:]
                    if (size(L,1)>1)
                        for ii=1:mesh.npoin
                            L[ip,ii] = 0.0
                        end
                        L[ip,ip] =1.0
                    end
                end
            end
        end
    else        
        for iedge = 1:size(mesh.bdy_edge_comp,1)
            iel = mesh.bdy_edge_in_elem[iedge]
            comp = mesh.bdy_edge_comp[iedge]
            if (mesh.bdy_edge_type[iedge] != "periodic1" && mesh.bdy_edge_type[iedge] !="periodic2")
                for k=1:mesh.ngl
                    tag = mesh.bdy_edge_type[iedge]
                    ip = mesh.poin_in_bdy_edge[iedge,k]
                    m=1
                    l=1
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            if (mesh.connijk[ii,jj,iel] == ip)
                                mm=jj
                                ll=ii
                            end
                        end
                    end
                    q[ip,:] .= 0.0
                    flux = zeros(size(q,2),1)
                    rhs[l,m,iel,:] .= rhs[l,m,iel,:] .+ flux[:]
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
end
