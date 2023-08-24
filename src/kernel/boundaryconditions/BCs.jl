#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
#include("../../io/mod_inputs.jl")
#include("../mesh/metric_terms.jl")
#include("../mesh/mesh.jl")
include("../operators/operators.jl")
#include("../abstractTypes.jl")
#include("../bases/basis_structs.jl")
#include("../infrastructure/element_matrices.jl")
include("../infrastructure/Kopriva_functions.jl")
#include("../infrastructure/2D_3D_structures.jl")
include("custom_bcs.jl")

function apply_periodicity!(SD::NSD_1D, rhs, qp, mesh, inputs, QT, metrics, ψ, dψ, ω, t, nvars)

    #NOTICE: " apply_periodicity!() in 1D is now only working for nvars=1!"
    
    #
    # 1D periodic
    #
    qp[mesh.npoin_linear] = 0.5*(qp[mesh.npoin_linear] .+ qp[1])
    qp[1] = qp[mesh.npoin_linear]

end

function apply_boundary_conditions!(SD::NSD_1D, rhs, qp, mesh,inputs, QT, metrics, ψ, dψ, ω, t, nvars;L=zeros(1,1))
    #If Neumann conditions are needed compute gradient
    calc_grad = false
    #   for key in keys(inputs)
    #     if (inputs[key] == "dirichlet" || inputs[key] == "neumann" || inputs[key] == "dirichlet/neumann")
    calc_grad = true
    #    end
    #  end
    gradq = zeros(mesh.npoin,nvars)
    #TODO remake build custom_bcs for new boundary data
    #if (calc_grad)
    #    gradq = build_gradient(SD, QT::Inexact, qp, ψ, dψ, ω, mesh, metrics,gradq,nvars)
    build_custom_bcs!(t,mesh,qp,gradq,rhs,SD,nvars,metrics,ω,dirichlet!,neumann,L,inputs)
    #end
end

#function apply_boundary_conditions!(SD::NSD_2D, rhs, qp, gradq, mesh, inputs, QT, metrics, ψ, dψ, ω, t, nvars;L=zeros(1,1))

function apply_boundary_conditions!(u, params, t)
    
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

    #gradq = zeros(2, 1, 1) #zeros(2,mesh.npoin,nvars)
    fill!(params.gradu, zero(Float64))
    
    #bdy_flux_q = zeros(mesh.ngl,nface,2,nvars)
    #exact = zeros(mesh.ngl,nface,nvars)
    #penalty =0.0#50000
    #nx = metrics.nx
    #ny = metrics.ny
    ##TODO remake build custom_bcs for new boundary data
    ##if (calc_grad)
    ##    gradq = build_gradient(SD, QT::Inexact, qp, ψ, dψ, ω, mesh, metrics,gradq,nvars)
    build_custom_bcs!(t, params.mesh, params.ubdy, params.uaux, params.gradu, params.bdy_flux,
                      params.rhs_el,
                      params.SD, params.neqs,
                      params.metrics, params.ω,
                      dirichlet!, neumann,
                      zeros(1,1), params.inputs)
    #end
    
end
#=

function build_custom_bcs!(t,mesh,q,gradq,rhs,::NSD_2D,nvars,metrics,ω,dirichlet!,neumann,L,inputs)

#error("QUI LinearCLaw/user_bc.jl")

c  = 1.0
x0 = y0 = -0.8
kx = ky = sqrt(2.0)/2
ω  = 0.2
d  = 0.5*ω/sqrt(log(2.0)); d2 = d*d
@inline f(x,y,t) = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)

#idx = (ieq-1)*mesh.npoin
#q[idx+1:ieq*mesh.npoin]
for ieq = 1:nvars
for iedge = 1:size(mesh.bdy_edge_comp,1)
iel = mesh.bdy_edge_in_elem[iedge]
comp = mesh.bdy_edge_comp[iedge]
for k=1:mesh.ngl
if (mesh.bdy_edge_type[iedge] != "periodic1" && mesh.bdy_edge_type[iedge] !="periodic2")
tag = mesh.bdy_edge_type[iedge]
ip = mesh.poin_in_bdy_edge[iedge,k]
idx = (ieq - 1)*mesh.npoin + ip

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
#if (inputs[:luser_bc])

e = f(x,y,t) #exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
q_bdyvalue = [e, kx*e/c, ky*e/c]                    
q[idx] = q_bdyvalue[ieq]

#q[ip,:] = dirichlet!(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag,inputs)
#bdy_flux = zeros(size(q,2),1)
#bdy_flux = (ω[k]*metrics.Jef[k,iedge]).*neumann(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag,inputs)
#END SM test

#q[ip,:] = dirichlet!(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag,inputs)
#bdy_flux = (ω[k]*metrics.Jef[k,iedge]).*neumann(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag,inputs)
#else
#    q[ip,:] .= 0.0
#    bdy_flux = zeros(size(q,2),1)
#end
#rhs[l,m,iel,:] .= rhs[l,m,iel,:] #.+ bdy_flux[:]
#  if (size(L,1)>1)
#      for ii=1:mesh.npoin
#          L[ip,ii] = 0.0
#      end
#      L[ip,ip] =1.0
#  end
end
end
end
end
end
=#

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
        qbdy = zeros(size(q,2),1) #NO
        
        qbdy[:] .= 4325789.0
        if (inputs[:luser_bc])
            qbdy = dirichlet!(q[ip,:], gradq[ip,:], x, t, mesh, metrics, "notag", qbdy, inputs)
            bdy_flux = (ω[k]*neumann(q[ip,:],gradq[ip,:],x,t,mesh,metrics,inputs))
        else
            q[ip,:] .= 0.0
            bdy_flux = zeros(size(q,2),1)
        end
        
        rhs[k,iel,:] .= rhs[k,iel,:] .+ bdy_flux[:]
        for var =1:size(q,2)
            if !(AlmostEqual(qbdy[var],4325789.0))
                #@info var,x,y,qbdy[var]
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



function build_custom_bcs!(t, mesh, qbdy, q, gradq, bdy_flux, rhs, ::NSD_2D, nvars, metrics, ω, dirichlet!, neumann, L, inputs)
    
    for iedge = 1:mesh.nedges_bdy 
        iel  = mesh.bdy_edge_in_elem[iedge]

        if mesh.bdy_edge_type[iedge] != "periodic1" && mesh.bdy_edge_type[iedge] != "periodic2"
            tag = mesh.bdy_edge_type[iedge]
            for k=1:mesh.ngl
                ip = mesh.poin_in_bdy_edge[iedge,k]
                
                fill!(qbdy, 4325789.0)
                ipp = 1 #ip               
                qbdy = dirichlet!(@view(qbdy[:]),@view(gradq[:,ipp,:]), mesh.x[ip], mesh.y[ip], t, tag, inputs)
                
                mm=1; ll=1
                for jj=1:mesh.ngl, ii=1:mesh.ngl
                    if (mesh.connijk[ii,jj,iel] == ip)
                        mm=jj
                        ll=ii
                    end
                end
                                
                for var =1:nvars
                    if !(AlmostEqual(qbdy[var],4325789.0)) # WHAT's this for?
                        q[ip,var]          = qbdy[var]
                        rhs[ll,mm,iel,var] = 0.0 #WHAT DOES THIS DO? here is only updated the  `ll` and `mm` row outside of any ll or mm loop
                    end
                end
            end
        end
    end
end

#=
function yt_build_custom_bcs!(t, mesh, qbdy, q, gradq, bdy_flux, rhs, ::NSD_2D, nvars, metrics, ω, dirichlet!, neumann, L, inputs)
    
    # nedges =  mesh.nedges_int
    # q_size = size(q, 2)
    # L_size = size(L, 1)
    # fill!(bdy_flux, zero(Float64))
        
    for iedge = 1:size(mesh.bdy_edge_in_elem,1)
        iel  = mesh.bdy_edge_in_elem[iedge]

        if mesh.bdy_edge_type[iedge] != "periodic1" && mesh.bdy_edge_type[iedge] != "periodic2"
            tag = mesh.bdy_edge_type[iedge]
            for k=1:mesh.ngl
                ip = mesh.poin_in_bdy_edge[iedge,k]
                ωJacedge = ω[k]*metrics.Jef[k,iedge]
                
                #bdy_flux = zeros(q_size,1)
                x    = mesh.x[ip]
                y    = mesh.y[ip]

                fill!(qbdy, 4325789.0)
                #flags = zeros(q_size,1) #NO, NEVER ALLOCATE INSIDE A LOOP
                if (inputs[:luser_bc]) #in contrast to what?
                    #q[ip,:], flags = dirichlet!(q[ip,:],gradq[:,ip,:],x,y,t,mesh,metrics,tag,qbdy)
                    ipp=1 #ip               
                    qbdy = dirichlet!(@view(q[ip,:]),@view(gradq[:,ipp,:]),x,y,t,mesh,metrics,tag,qbdy,inputs)
                    ##SM change this to set bdy_flux to zero and do not allocate gradq unless neumann is required explicitly by the user
                    bdy_flux .= ωJacedge.*neumann(q[ip,:],gradq[:,ipp,:],x,y,t,mesh,metrics,tag,inputs)
                else
                    q[ip,:] .= 0.0
                end

                
                mm=1; ll=1
                for jj=1:mesh.ngl
                    for ii=1:mesh.ngl
                        if (mesh.connijk[ii,jj,iel] == ip)
                            mm=jj
                            ll=ii
                        end
                    end
                end
                rhs[ll,mm,iel,:] .= rhs[ll,mm,iel,:] .+ bdy_flux[:]
                
                
                for var =1:nvars
                    if !(AlmostEqual(qbdy[var],4325789.0)) # WHAT's this for?
                        #@info var,x,y,qbdy[var]
                        
                        rhs[ll,mm,iel,var] = 0.0 #WHAT DOES THIS DO? here is only updated the  `ll` and `mm` row outside of any ll or mm loop
                        
                        q[ip,var]          = qbdy[var]
                    end
                end
                # if (L_size > 1)  ALL THI SABOUT L should exist when we use L
                #, which we don't ever unless we explicitly build it.
                #     for ii=1:mesh.npoin
                #         L[ip,ii] = 0.0
                #     end
                #     L[ip,ip] = 1.0
                # end
                
            end
        end
    end
end
=#
