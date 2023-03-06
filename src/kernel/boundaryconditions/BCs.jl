#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("../../io/mod_inputs.jl")
include("../operators/operators.jl")
include("../AbstractTypes.jl")
include("../bases/basis_structs.jl")
include("../infrastructure/element_matrices.jl")
include("../infrastructure/Kopriva_functions.jl")
include("../infrastructure/2D_3D_structures.jl")
include("../mesh/metric_terms.jl")
include("../mesh/mesh.jl")


function apply_periodicity!(rhs,qp,mesh,inputs, SD::NSD_1D,QT,metrics,ψ,dψ, ω,t,BCT,nvars)
    
    if (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="periodic" || haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="periodic")
        #
        # B.C.: 1D periodic
        #
        qp[mesh.npoin_linear,:] .= 0.5*(qp[mesh.npoin_linear,:] .+ qp[1,:])
        qp[1,:] .= qp[mesh.npoin_linear,:]
    elseif (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="dirichlet" || haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="dirichlet")
        #
        # B.C.: solid wall
        #
        qp[1] = 0.0
        qp[mesh.npoin_linear] = 0.0
    end
end

function apply_boundary_conditions!(rhs,qp,mesh,inputs, SD::NSD_1D,QT,metrics,ψ,dψ, ω,t,BCT,nvars)
    nothing
end
    
function apply_periodicity!(rhs,qp,mesh,inputs, SD::NSD_2D,QT,metrics,ψ,dψ, ω,t,BCT,nvars)
   #Periodic boundaries
   if (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="periodic" || haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="periodic")
      inputs[:xmax_bc] = "periodic"
      inputs[:xmin_bc] = "periodic"
      # B.C.: periodic
      for I in keys(mesh.xperiodicity)
         ip = mesh.xperiodicity[I]
         qp[I,:] .= 0.5*(qp[I,:] .+ qp[ip,:])
         qp[ip,:] .= qp[I,:]
      end
    end 
    if (haskey(inputs, :ymin_bc) && inputs[:ymin_bc]=="periodic") || (haskey(inputs, :ymax_bc) && inputs[:ymax_bc]=="periodic")
        inputs[:ymax_bc] = "periodic"
        inputs[:ymin_bc] = "periodic"
        for I in keys(mesh.yperiodicity)
            ip = mesh.yperiodicity[I]
            qp[I,:] .= 0.5*(qp[I,:] .+ qp[ip,:])
            qp[ip,:] .= qp[I,:]
        end
    end
end
function apply_boundary_conditions!(rhs,qp,mesh,inputs, SD::NSD_2D,QT,metrics,ψ,dψ, ω,t,BCT,nvars)
   #If Neumann conditions are needed compute gradient
   calc_grad = false
   for key in keys(inputs)
      if (inputs[key] == "dirichlet" || inputs[key] == "neumann" || inputs[key] == "dirichlet/neumann")
          calc_grad = true
      end
   end
   nface = size(metrics.Jef,2)
   dqdx_st = zeros(nvars,2)
   q_st = zeros(nvars,1)
   gradq = zeros(2,mesh.npoin,nvars)
   flux_q = zeros(mesh.ngl,nface,2,nvars)
   exact = zeros(mesh.ngl,nface,nvars)
   penalty =160.0#50000
   nx = zeros(mesh.ngl,nface)
   ny = zeros(mesh.ngl,nface)
   if (calc_grad)
       gradq = build_gradient(SD, QT::Inexact, qp, ψ, dψ, ω, mesh, metrics,gradq,nvars)
       build_custom_bcs!(t,mesh,qp,gradq,SD,BCT,exact,flux_q,nx,ny,nvars)
   end
   #Dirichlet/Neumann boundaries using SIPG
   # NOTE We do not need to compute a RHS contribution for the Right element as it represents the outside of the computational domain here we only compute it's effect on the Left element
   if (haskey(inputs, :xmin_bc) && (inputs[:xmin_bc]!="periodic"))
#      @info "applying sipg" 
      for iface=1:size(mesh.xmin_faces,2)
          for k=1:mesh.ngl
              for i=1:mesh.ngl
                  iel = mesh.xmin_facetoelem[iface] 
                  mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,iface]/metrics.Je[1,k,iel]/2
                  ip = mesh.xmin_faces[k,iface]
                  #@info "1",rhs[1,k,iel,1:nvars], mesh.x[ip],mesh.y[ip]
                  dqdx_st[:,1] .= 0.5*(gradq[1,ip,:] .+ flux_q[k,iface,1,:] .- nx[k,iface]*mu.*(exact[k,iface,:].-qp[ip,1:nvars]))
                  dqdx_st[:,2] .= 0.5*(gradq[2,ip,:] .+ flux_q[k,iface,2,:] .- ny[k,iface]*mu.*(exact[k,iface,:].-qp[ip,1:nvars]))
                  q_st[:] .= 0.5*(qp[ip,1:nvars] + exact[k,iface,:]) 
                  rhs[1,k,iel,1:nvars] .-= ω[k]*metrics.Jef[k,iface]*ψ[i,k].*(nx[k,iface]*dqdx_st[1:nvars,1]) #.+ ny[ip]*dqdx_st[1:nvars,2])
                  #@info "2",rhs[1,k,iel,1:nvars],mesh.x[ip],mesh.y[ip]
                  rhs[1,k,iel,1:nvars] .-= ω[k]*metrics.Jef[k,iface]*(nx[k,iface])*dψ[i,k].*(qp[ip,1:nvars] .- q_st[1:nvars])
                  #@info "3",rhs[1,k,iel,1:nvars],mesh.x[ip],mesh.y[ip]
                  #end
              end
          end
      end
  end
          
  if (haskey(inputs, :xmax_bc) && (inputs[:xmax_bc]!="periodic"))
            #Right boundary
 #@info "applying sipg"
      disp = size(mesh.xmin_faces,2)
      for iface=1:size(mesh.xmax_faces,2)
          for k=1:mesh.ngl
              for i=1:mesh.ngl
                  ip = mesh.xmax_faces[k,iface]
                  iel = mesh.xmax_facetoelem[iface] 
                  mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,iface+size(mesh.xmin_faces,2)]/metrics.Je[mesh.ngl,k,iel]/2
                  dqdx_st[:,1] .= 0.5*(gradq[1,ip,:] .+ flux_q[k,iface+disp,1,:] - nx[k,iface+disp]*mu.*(exact[k,iface+disp,:].-qp[ip,1:nvars]))
                  dqdx_st[:,2] .= 0.5*(gradq[2,ip,:] .+ flux_q[k,iface+disp,2,:] - ny[k,iface+disp]*mu.*(exact[k,iface+disp,:].-qp[ip,1:nvars]))
                  q_st[:] = 0.5*(qp[ip,1:nvars] + exact[k,iface+disp,:])
                  rhs[mesh.ngl,k,iel,1:nvars] .-= ω[k]*metrics.Jef[k,iface+size(mesh.xmin_faces,2)]*ψ[i,k].*(nx[k,iface+disp]*dqdx_st[1:nvars,1]) #.+ ny[ip]*dqdx_st[1:nvars,2])
                  rhs[mesh.ngl,k,iel,1:nvars] .-= ω[k]*metrics.Jef[k,iface+size(mesh.xmin_faces,2)]*(nx[k,iface+disp])*dψ[i,k].*(qp[ip,1:nvars] .- q_st[1:nvars])
                  #if (qp[ip,1] < -0.01)
#                      @info rhs[ip],ω[k],metrics.Jef[k,iface],(nx[ip]),dψ[i,k],(qp[ip,1] - q_st[1]),q_st[1]
                  #end 
              end
          end
      end
  end
  if (haskey(inputs, :ymin_bc) && (inputs[:ymin_bc]!="periodic"))
         #bottom boundary
       disp = size(mesh.xmin_faces,2)+size(mesh.xmax_faces,2)
       for iface=1:size(mesh.ymin_faces,2)
          for k=1:mesh.ngl
              for i=1:mesh.ngl
                  ip = mesh.ymin_faces[k,iface]
                  iel = mesh.ymin_facetoelem[iface] 
                  mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,iface+disp]/metrics.Je[k,1,iel]/2
                  dqdx_st[:,1] .= 0.5*(gradq[1,ip,:] .+ flux_q[k,iface+disp,1,:] - nx[k,iface+disp]*mu.*(exact[k,iface+disp,:].-qp[ip,1:nvars]))
                  dqdx_st[:,2] .= 0.5*(gradq[2,ip,:] .+ flux_q[k,iface+disp,2,:] - ny[k,iface+disp]*mu.*(exact[k,iface+disp,:].-qp[ip,1:nvars]))
                  q_st[:] = 0.5*(qp[ip,1:nvars] + exact[k,iface+disp,:])
                #  @info "1",rhs[k,1,iel,1:nvars], mesh.x[ip],mesh.y[ip]
                   rhs[k,1,iel,1:nvars] .-= ω[k]*metrics.Jef[k,iface+disp]*ψ[i,k]*(ny[k,iface+disp]*dqdx_st[1:nvars,2])
                #  @info "2",rhs[k,1,iel,1:nvars], mesh.x[ip],mesh.y[ip]
                  rhs[k,1,iel,1:nvars] .-= ω[k]*metrics.Jef[k,iface+disp]*(ny[k,iface+disp])*dψ[i,k].*(qp[ip,1:nvars] .- q_st[1:nvars])
                #  @info "3",rhs[k,1,iel,1:nvars], mesh.x[ip],mesh.y[ip]
                   #if (qp[ip,1] < -0.01)
              #   @info "ymin", rhs[k,1,iel,3],rhs[k,1,iel,2]
#     @info "ymin",mesh.x[mesh.connijk[k,1,iel]],mesh.y[mesh.connijk[k,1,iel]],rhs[k,1,iel,3],ψ[i,k]*(ny[k,iface+disp]*dqdx_st[3,2]),(ny[k,iface+disp])*dψ[i,k]*(qp[ip,3] - q_st[3]),q_st[3]
                  #end 
              end
          end
      end

   end
   if (haskey(inputs, :ymax_bc) && (inputs[:ymax_bc]!="periodic"))
            #top boundary
      disp = size(mesh.xmin_faces,2)+size(mesh.xmax_faces,2) + size(mesh.ymin_faces,2)
      for iface=1:size(mesh.ymax_faces,2)
          for k=1:mesh.ngl
              for i=1:mesh.ngl
                   ip = mesh.ymax_faces[k,iface]
                   iel = mesh.ymax_facetoelem[iface] 
                   mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,iface+disp]/metrics.Je[k,mesh.ngl,iel]/2
                   dqdx_st[:,1] .= 0.5*(gradq[1,ip,:] .+ flux_q[k,iface+disp,1,:] - nx[k,iface+disp]*mu.*(exact[k,iface+disp,:].-qp[ip,1:nvars]))
                   dqdx_st[:,2] .= 0.5*(gradq[2,ip,:] .+ flux_q[k,iface+disp,2,:] - ny[k,iface+disp]*mu.*(exact[k,iface+disp,:].-qp[ip,1:nvars]))
                   q_st[:] = 0.5*(qp[ip,1:nvars] + exact[k,iface+disp,:])
                   rhs[k,mesh.ngl,iel,1:nvars] .-= ω[k]*metrics.Jef[k,iface+disp]*ψ[i,k].*(0.0*dqdx_st[1:nvars,1] .+ ny[k,iface+disp]*dqdx_st[1:nvars,2])
                   rhs[k,mesh.ngl,iel,1:nvars] .-= ω[k]*metrics.Jef[k,iface+disp]*(ny[k,iface+disp])*dψ[i,k].*(qp[ip,1:nvars] .- q_st[1:nvars])
                  #if (qp[ip,1] < -0.01)
                   #   @info "ymax",nx[ip],ny[ip],rhs[k,mesh.ngl,iel],ω[k]*metrics.Jef[k,iface]*ψ[i,k]*(ny[ip]*dqdx_st[1,2]),ω[k]*metrics.Jef[k,iface]*(ny[ip])*dψ[i,k]*(qp[ip,1] - q_st[1])
                  #end  
             end
          end
      end
   end
end
function apply_periodicity!(qp,mesh,inputs,SD::NSD_3D,QT,metrics,ψ,dψ, ω,t,BCT,nvars)

    if (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="periodic" || haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="periodic")
        inputs[:xmax_bc] = "periodic"
        inputs[:xmin_bc] = "periodic"
        # B.C.: periodic
        for I in keys(mesh.xperiodicity)
            ip = mesh.xperiodicity[I]
            qp[I,1:nvars] .= 0.5*(qp[I,1:nvars] .+ qp[ip,1:nvars])
            qp[ip,1:nvars] .= qp[I,1:nvars]
        end
    end
    if (haskey(inputs, :ymin_bc) && inputs[:ymin_bc]=="periodic") || (haskey(inputs, :ymax_bc) && inputs[:ymax_bc]=="periodic")
        inputs[:ymax_bc] = "periodic"
        inputs[:ymin_bc] = "periodic"
        for I in keys(mesh.yperiodicity)
            ip = mesh.yperiodicity[I]
            qp[I,:] .= 0.5*(qp[I,1:nvars] .+ qp[ip,1:nvars])
            qp[ip,:] .= qp[I,:]
        end
    end
    if (haskey(inputs, :zmin_bc) && inputs[:zmin_bc]=="periodic") || (haskey(inputs, :zmax_bc) && inputs[:zmax_bc]=="periodic")
        inputs[:zmax_bc] = "periodic"
        inputs[:zmin_bc] = "periodic"
        for I in keys(mesh.zperiodicity)
            ip = mesh.zperiodicity[I]
            qp[I,:] .= 0.5*(qp[I,1:nvars] .+ qp[ip,1:nvars])
            qp[ip,:] .= qp[I,:]
        end
    end
end
function apply_boundary_conditions!(rhs,qp,mesh,inputs, SD::NSD_3D,QT,metrics,ψ,dψ, ω,t,BCT,nvars)
   calc_grad = false
   for key in keys(inputs)
      if (inputs[key] == "dirichlet" || inputs[key] == "neumann" || inputs[key] == "dirichlet/neumann")
          calc_grad = true
      end
   end
#  @info calc_grad
   dqdx_st = zeros(nvars,3)
   q_st = zeros(nvars,1)
   gradq = zeros(mesh.npoin,3,nvars)
   flux_q = zeros(mesh.npoin,3,nvars)
   exact = zeros(mesh.npoin,nvars)
   penalty = 0.0#50000
   nx = zeros(mesh.npoin,1)
   ny = zeros(mesh.npoin,1)
   nz = zeors(mesh.npoin,1)
   if (calc_grad)
       gradq = build_gradient(SD, QT::Inexact, qp, ψ, dψ, ω, mesh, metrics)
       build_custom_bcs!(t,mesh,qp,gradq,SD,BCT,exact,flux_q,nx,ny)
   end
   #Dirichlet/Neumann boundaries using SIPG
   # NOTE We do not need to compute a RHS contribution for the Right element as it represents the outside of the computational domain here we only compute it's effect on the Left element
   if (haskey(inputs, :xmin_bc) && (inputs[:xmin_bc]!="periodic"))
#      @info "applying sipg"
      for iface=1:size(mesh.xmin_faces,2)
          for k=1:mesh.ngl
              for l=1:mesh.ngl
                  for i=1:mesh.ngl
                      for j=1:mesh.ngl
                          iel = mesh.xmin_facetoelem[iface]
                          mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,l,iel]/metrics.Je[1,k,l,iface]/2
                          ip = mesh.xmin_faces[k,l,iface]
                          dqdx_st[:,1] = 0.5*(gradq[1,ip,:] .+ flux_q[ip,1,:] - nx[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,2] = 0.5*(gradq[2,ip,:] .+ flux_q[ip,2,:] - ny[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,3] = 0.5*(gradq[3,ip,:] .+ flux_q[ip,3,:] - nz[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          q_st[:] = 0.5*(qp[ip,1:nvars] + exact[ip,:])
                          rhs[1,i,j,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*ψ[i,k]*ψ[j,l]*(nx[ip]*dqdx_st[1:nvars,1] .+ ny[ip]*dqdx_st[1:nvars,3] .+ nz[ip]*dqdx_st[1:nvars,3])
                          rhs[1,i,j,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*(nx[ip])*dψ[i,k]*dψ[j,l]*(qp[ip,1:nvars] - q_st[1:nvars])
                  #if (qp[ip,1] < -0.01)
                   # @info calc_grad,rhs[ip],ω[k],metrics.Jef[k,iface],(nx[ip]),dψ[i,k],(qp[ip,1] - q_st[1])
                  #end
                      end
                  end
              end
          end
      end
  end
  if (haskey(inputs, :xmax_bc) && (inputs[:xmax_bc]!="periodic"))
#      @info "applying sipg"
      disp = size(mesh.xmin_faces,2)
      for iface=1:size(mesh.xmax_faces,2)
          for k=1:mesh.ngl
              for l=1:mesh.ngl
                  for i=1:mesh.ngl
                      for j=1:mesh.ngl
                          iel = mesh.xmax_facetoelem[iface]
                          mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,l,iel]/metrics.Je[mesh.ngl,k,l,iface]/2
                          ip = mesh.xmax_faces[k,l,iface]
                          dqdx_st[:,1] = 0.5*(gradq[1,ip,:] .+ flux_q[ip,1,:] - nx[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,2] = 0.5*(gradq[2,ip,:] .+ flux_q[ip,2,:] - ny[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,3] = 0.5*(gradq[3,ip,:] .+ flux_q[ip,3,:] - nz[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          q_st[:] = 0.5*(qp[ip,1:nvars] + exact[ip,:])
                          rhs[mesh.ngl,i,j,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*ψ[i,k]*ψ[j,l]*(nx[ip]*dqdx_st[1:nvars,1] .+ ny[ip]*dqdx_st[1:nvars,3] .+ nz[ip]*dqdx_st[1:nvars,3])
                          rhs[mesh.ngl,i,j,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*(nx[ip])*dψ[i,k]*dψ[j,l]*(qp[ip,1:nvars] - q_st[1:nvars])      
           #if (qp[ip,1] < -0.01)
                   # @info calc_grad,rhs[ip],ω[k],metrics.Jef[k,iface],(nx[ip]),dψ[i,k],(qp[ip,1] - q_st[1])
                  #end
                      end
                  end
              end
          end
      end
  end
  if (haskey(inputs, :ymin_bc) && (inputs[:ymin_bc]!="periodic"))
#      @info "applying sipg"
      disp = disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2)
      for iface=1:size(mesh.ymin_faces,2)
          for k=1:mesh.ngl
              for l=1:mesh.ngl
                  for i=1:mesh.ngl
                      for j=1:mesh.ngl
                          iel = mesh.ymin_facetoelem[iface]
                          mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,l,iel]/metrics.Je[k,1,l,iface]/2
                          ip = mesh.ymin_faces[k,l,iface]
                          dqdx_st[:,1] = 0.5*(gradq[1,ip,:] .+ flux_q[ip,1,:] - nx[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,2] = 0.5*(gradq[2,ip,:] .+ flux_q[ip,2,:] - ny[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,3] = 0.5*(gradq[3,ip,:] .+ flux_q[ip,3,:] - nz[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          q_st[:] = 0.5*(qp[ip,1:nvars] + exact[ip,:])
                          rhs[i,1,j,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*ψ[i,k]*ψ[j,l]*(nx[ip]*dqdx_st[1:nvars,1] .+ ny[ip]*dqdx_st[1:nvars,3] .+ nz[ip]*dqdx_st[1:nvars,3])
                          rhs[i,1,j,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*(nx[ip])*dψ[i,k]*dψ[j,l]*(qp[ip,1:nvars] - q_st[1:nvars])
                   #if (qp[ip,1] < -0.01)
                   # @info calc_grad,rhs[ip],ω[k],metrics.Jef[k,iface],(nx[ip]),dψ[i,k],(qp[ip,1] - q_st[1])
                  #end
                      end
                  end
              end
          end
      end
  end
  if (haskey(inputs, :ymax_bc) && (inputs[:ymax_bc]!="periodic"))
#      @info "applying sipg"
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2) + size(mesh.ymin_faces,2)
      for iface=1:size(mesh.ymax_faces,2)
          for k=1:mesh.ngl
              for l=1:mesh.ngl
                  for i=1:mesh.ngl
                      for j=1:mesh.ngl
                          iel = mesh.ymax_facetoelem[iface]
                          mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,l,iel]/metrics.Je[k,mesh.ngl,l,iface]/2
                          ip = mesh.ymax_faces[k,l,iface]
                          dqdx_st[:,1] = 0.5*(gradq[1,ip,:] .+ flux_q[ip,1,:] - nx[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,2] = 0.5*(gradq[2,ip,:] .+ flux_q[ip,2,:] - ny[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,3] = 0.5*(gradq[3,ip,:] .+ flux_q[ip,3,:] - nz[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          q_st[:] = 0.5*(qp[ip,1:nvars] + exact[ip,:])
                          rhs[i,mesh.ngl,j,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*ψ[i,k]*ψ[j,l]*(nx[ip]*dqdx_st[1:nvars,1] .+ ny[ip]*dqdx_st[1:nvars,3] .+ nz[ip]*dqdx_st[1:nvars,3])
                          rhs[i,mesh.ngl,j,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*(nx[ip])*dψ[i,k]*dψ[j,l]*(qp[ip,1:nvars] - q_st[1:nvars])
                  #if (qp[ip,1] < -0.01)
                   # @info calc_grad,rhs[ip],ω[k],metrics.Jef[k,iface],(nx[ip]),dψ[i,k],(qp[ip,1] - q_st[1])
                  #end
                      end
                  end
              end
          end
      end
  end
  if (haskey(inputs, :zmin_bc) && (inputs[:zmin_bc]!="periodic"))
#      @info "applying sipg"
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2) + size(mesh.ymin_faces,2) + size(mesh.ymax_faces,2)
      for iface=1:size(mesh.zmin_faces,2)
          for k=1:mesh.ngl
              for l=1:mesh.ngl
                  for i=1:mesh.ngl
                      for j=1:mesh.ngl
                          iel = mesh.zmin_facetoelem[iface]
                          mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,l,iel]/metrics.Je[k,l,1,iface]/2
                          ip = mesh.zmin_faces[k,l,iface]
                          dqdx_st[:,1] = 0.5*(gradq[1,ip,:] .+ flux_q[ip,1,:] - nx[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,2] = 0.5*(gradq[2,ip,:] .+ flux_q[ip,2,:] - ny[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,3] = 0.5*(gradq[3,ip,:] .+ flux_q[ip,3,:] - nz[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          q_st[:] = 0.5*(qp[ip,1:nvars] + exact[ip,:])
                          rhs[i,j,1,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*ψ[i,k]*ψ[j,l]*(nx[ip]*dqdx_st[1:nvars,1] .+ ny[ip]*dqdx_st[1:nvars,3] .+ nz[ip]*dqdx_st[1:nvars,3])
                          rhs[i,j,1,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*(nx[ip])*dψ[i,k]*dψ[j,l]*(qp[ip,1:nvars] - q_st[1:nvars])
                  #if (qp[ip,1] < -0.01)
                   # @info calc_grad,rhs[ip],ω[k],metrics.Jef[k,iface],(nx[ip]),dψ[i,k],(qp[ip,1] - q_st[1])
                  #end
                      end
                  end 
              end
          end
      end
  end
  if (haskey(inputs, :zmax_bc) && (inputs[:zmax_bc]!="periodic"))
#      @info "applying sipg"
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2) + size(mesh.ymin_faces,2) + size(mesh.ymax_faces,2) + size(mesh.zmin_faces,2)
      for iface=1:size(mesh.zmax_faces,2)
          for k=1:mesh.ngl
              for l=1:mesh.ngl
                  for i=1:mesh.ngl
                      for j=1:mesh.ngl
                          iel = mesh.zmax_facetoelem[iface]
                          mu = penalty * (mesh.ngl)*(mesh.ngl-1)*metrics.Jef[k,l,iel]/metrics.Je[k,l,mesh.ngl,iface]/2
                          ip = mesh.zmax_faces[k,l,iface]
                          dqdx_st[:,1] = 0.5*(gradq[1,ip,:] .+ flux_q[ip,1,:] - nx[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,2] = 0.5*(gradq[2,ip,:] .+ flux_q[ip,2,:] - ny[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          dqdx_st[:,3] = 0.5*(gradq[3,ip,:] .+ flux_q[ip,3,:] - nz[ip]*mu.*(exact[ip,:].-qp[ip,1:nvars]))
                          q_st[:] = 0.5*(qp[ip,1:nvars] + exact[ip,:])
                          rhs[i,j,mesh.ngl,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*ψ[i,k]*ψ[j,l]*(nx[ip]*dqdx_st[1:nvars,1] .+ ny[ip]*dqdx_st[1:nvars,3] .+ nz[ip]*dqdx_st[1:nvars,3])
                          rhs[i,j,mesh.ngl,iel,1:nvars] .-= ω[k]*ω[l]*metrics.Jef[k,l,iface]*(nx[ip])*dψ[i,k]*dψ[j,l]*(qp[ip,1:nvars] - q_st[1:nvars])
                  #if (qp[ip,1] < -0.01)
                   # @info calc_grad,rhs[ip],ω[k],metrics.Jef[k,iface],(nx[ip]),dψ[i,k],(qp[ip,1] - q_st[1])
                  #end
                      end
                  end
              end
          end
      end
  end
end

function build_custom_bcs!(t,mesh,q,gradq,::NSD_2D,::DefaultBC,exact_q,flux_q,nx,ny,nvars)
    for k=1:size(mesh.xmin_faces,2)
          for i=1:mesh.ngl
              nx[i,k] =1.0
              ny[i,k] =0.0
              ip = mesh.xmin_faces[i,k]
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,k,1] = -q[ip,1]
              exact_q[i,k,2] = q[ip,2] - 2*unl*nx[i,k]
              exact_q[i,k,3] = q[ip,3] - 2*unl*ny[i,k]
              for var=1:nvars
                  unl = nx[i,k]*gradq[1,ip,var]+ny[i,k]*gradq[2,ip,var]
                  flux_q[i,k,1,var] = gradq[1,ip,var] - 2*gradq[1,ip,var]*unl
                  flux_q[i,k,2,var] = gradq[2,ip,var] - 2*gradq[2,ip,var]*unl
              end
          end
      end
      disp = size(mesh.xmin_faces,2)
      for k=1:size(mesh.xmax_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =-1.0
              ny[i,kf] =0.0

              ip = mesh.xmax_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              exact_q[i,kf,1] = -q[ip,1]
              exact_q[i,kf,2] = q[ip,2] - 2*unl*nx[i,kf]
              exact_q[i,kf,3] = q[ip,3] - 2*unl*ny[i,kf]
              for var=1:nvars
                  unl = nx[i,kf]*gradq[1,ip,var]+ny[i,kf]*gradq[2,ip,var]
                  flux_q[i,kf,1,var] = gradq[1,ip,var] - 2*gradq[1,ip,var]*unl
                  flux_q[i,kf,2,var] = gradq[2,ip,var] - 2*gradq[2,ip,var]*unl
              end
          end
      end  
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2)
      for k=1:size(mesh.ymin_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =0.0
              ny[i,kf] =1.0

              ip = mesh.ymin_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              exact_q[i,kf,1] = -q[ip,1]
              exact_q[i,kf,2] = q[ip,2] - 2*unl*nx[i,kf]
              exact_q[i,kf,3] = q[ip,3] - 2*unl*ny[i,kf]
              for var=1:nvars
                  unl = nx[i,kf]*gradq[1,ip,var]+ny[i,kf]*gradq[2,ip,var]
                  flux_q[i,kf,1,var] = gradq[1,ip,var] - 2*gradq[1,ip,var]*unl
                  flux_q[i,kf,2,var] = gradq[2,ip,var] - 2*gradq[2,ip,var]*unl
              end
          end
      end
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2) + size(mesh.ymin_faces,2)
      for k=1:size(mesh.ymax_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =0.0
              ny[i,kf] =-1.0

              ip = mesh.ymax_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              exact_q[i,kf,1] = -q[ip,1]
              exact_q[i,kf,2] = q[ip,2] - 2*unl*nx[i,kf]
              exact_q[i,kf,3] = q[ip,3] - 2*unl*ny[i,kf]
              for var=1:nvars
                  unl = nx[i,kf]*gradq[1,ip,var]+ny[i,kf]*gradq[2,ip,var]
                  flux_q[i,kf,1,var] = gradq[1,ip,var] - 2*gradq[1,ip,var]*unl
                  flux_q[i,kf,2,var] = gradq[2,ip,var] - 2*gradq[2,ip,var]*unl
              end
          end
      end 
end

function build_custom_bcs!(t,mesh,q,gradq,::NSD_2D,::LinearClaw_KopNR,exact_q,flux_q,nx,ny,nvars)
      c  = 1.0
    x0 = y0 = -0.8
    kx = ky = sqrt(2)/2
    ω  = 0.2
    d  = 0.5*ω/sqrt(log(2)); d2 = d*d
      for k=1:size(mesh.xmin_faces,2)
          for i=1:mesh.ngl
              nx[i,k] =1.0
              ny[i,k] =0.0
              ip = mesh.xmin_faces[i,k]
              x = mesh.x[ip]
              y = mesh.y[ip]
              e = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,k,1] =  e
              exact_q[i,k,2] = kx*e/c
              exact_q[i,k,3] = ky*e/c
              flux_q[i,k,1,1] = -((2*kx^2)*x+(2*kx*ky*(y-y0)-(2*kx^2)*x0-2*kx*c*t))*e
              flux_q[i,k,2,1] = -((2*ky^2)*y+(2*kx*ky*(x-x0)-(2*ky^2)*y0-2*ky*c*t))*e
              flux_q[i,k,1,2] = flux_q[i,k,1,1]*kx/c
              flux_q[i,k,2,2] = flux_q[i,k,2,1]*kx/c
              flux_q[i,k,1,3] = flux_q[i,k,1,1]*ky/c
              flux_q[i,k,2,3] = flux_q[i,k,2,1]*ky/c
          end
      end
      disp = size(mesh.xmin_faces,2)
      for k=1:size(mesh.xmax_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =-1.0
              ny[i,kf] =0.0
             
              ip = mesh.xmax_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              x = mesh.x[ip]
              y = mesh.y[ip]
              e = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,kf,1] = e
              exact_q[i,kf,2] = kx*e/c
              exact_q[i,kf,3] = ky*e/c
              flux_q[i,kf,1,1] = -((2*kx^2)*x+(2*kx*ky*(y-y0)-(2*kx^2)*x0-2*kx*c*t))*e
              flux_q[i,kf,2,1] = -((2*ky^2)*y+(2*kx*ky*(x-x0)-(2*ky^2)*y0-2*ky*c*t))*e
              flux_q[i,kf,1,2] = flux_q[i,kf,1,1]*kx/c
              flux_q[i,kf,2,2] = flux_q[i,kf,2,1]*kx/c
              flux_q[i,kf,1,3] = flux_q[i,kf,1,1]*ky/c
              flux_q[i,kf,2,3] = flux_q[i,kf,2,1]*ky/c
          end
      end
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2)
      for k=1:size(mesh.ymin_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =0.0
              ny[i,kf] =-1.0
             
              ip = mesh.ymin_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              x = mesh.x[ip]
              y = mesh.y[ip]
              e = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,kf,1] = e
              exact_q[i,kf,2] = kx*e/c
              exact_q[i,kf,3] = ky*e/c
              flux_q[i,kf,1,1] = -((2*kx^2)*x+(2*kx*ky*(y-y0)-(2*kx^2)*x0-2*kx*c*t))*e
              flux_q[i,kf,2,1] = -((2*ky^2)*y+(2*kx*ky*(x-x0)-(2*ky^2)*y0-2*ky*c*t))*e
              flux_q[i,kf,1,2] = flux_q[i,kf,1,1]*kx/c
              flux_q[i,kf,2,2] = flux_q[i,kf,2,1]*kx/c
              flux_q[i,kf,1,3] = flux_q[i,kf,1,1]*ky/c
              flux_q[i,kf,2,3] = flux_q[i,kf,2,1]*ky/c
          end
      end
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2) + size(mesh.ymin_faces,2)
      for k=1:size(mesh.ymax_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =0.0
              ny[i,kf] =1.0

              ip = mesh.ymax_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              x = mesh.x[ip]
              y = mesh.y[ip]
              e = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,kf,1] = e
              exact_q[i,kf,2] = kx*e/c
              exact_q[i,kf,3] = ky*e/c
              flux_q[i,kf,1,1] = -((2*kx^2)*x+(2*kx*ky*(y-y0)-(2*kx^2)*x0-2*kx*c*t))*e 
              flux_q[i,kf,2,1] = -((2*ky^2)*y+(2*kx*ky*(x-x0)-(2*ky^2)*y0-2*ky*c*t))*e
              flux_q[i,kf,1,2] = flux_q[i,kf,1,1]*kx/c
              flux_q[i,kf,2,2] = flux_q[i,kf,2,1]*kx/c
              flux_q[i,kf,1,3] = flux_q[i,kf,1,1]*ky/c
              flux_q[i,kf,2,3] = flux_q[i,kf,2,1]*ky/c 
          end
      end
  #nx .= -nx
  #ny .= -ny
end

function build_custom_bcs!(t,mesh,q,gradq,::NSD_2D,::LinearClaw_KopRefxmax,exact_q,flux_q,nx,ny,nvars)
      c  = 1.0
    x0 = y0 = -0.8
    kx = ky = sqrt(2)/2
    ω  = 0.2
    d  = 0.5*ω/sqrt(log(2)); d2 = d*d
      for k=1:size(mesh.xmin_faces,2)
          for i=1:mesh.ngl
              nx[i,k] =1.0
              ny[i,k] =0.0
              ip = mesh.xmin_faces[i,k]
              x = mesh.x[ip]
              y = mesh.y[ip]
              e = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,k,1] =  e
              exact_q[i,k,2] = kx*e/c
              exact_q[i,k,3] = ky*e/c
              flux_q[i,k,1,1] = -((2*kx^2)*x+(2*kx*ky*(y-y0)-(2*kx^2)*x0-2*kx*c*t))*e
              flux_q[i,k,2,1] = -((2*ky^2)*y+(2*kx*ky*(x-x0)-(2*ky^2)*y0-2*ky*c*t))*e
              flux_q[i,k,1,2] = flux_q[i,k,1,1]*kx/c
              flux_q[i,k,2,2] = flux_q[i,k,2,1]*kx/c
              flux_q[i,k,1,3] = flux_q[i,k,1,1]*ky/c
              flux_q[i,k,2,3] = flux_q[i,k,2,1]*ky/c
          end
      end
      disp = size(mesh.xmin_faces,2)
      for k=1:size(mesh.xmax_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =-1.0
              ny[i,kf] =0.0

              ip = mesh.xmax_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              x = mesh.x[ip]
              y = mesh.y[ip]
              e = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,kf,1] = q[ip,1]
              exact_q[i,kf,2] = nx[i,kf]*kx*e/c
              exact_q[i,kf,3] = -nx[i,kf]*ky*e/c
              flux_q[i,kf,1,1] = gradq[1,ip,1]#-((2*kx^2)*x+(2*kx*ky*(y-y0)-(2*kx^2)*x0-2*kx*c*t))*e
              flux_q[i,kf,2,1] = gradq[1,ip,2]#-((2*ky^2)*y+(2*kx*ky*(x-x0)-(2*ky^2)*y0-2*ky*c*t))*e
              flux_q[i,kf,1,2] = nx[i,kf]*flux_q[i,kf,1,1]*kx/c
              flux_q[i,kf,2,2] = nx[i,kf]*flux_q[i,kf,2,1]*kx/c
              flux_q[i,kf,1,3] = -nx[i,kf]*flux_q[i,kf,1,1]*ky/c
              flux_q[i,kf,2,3] = -nx[i,kf]*flux_q[i,kf,2,1]*ky/c
          end
      end
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2)
      for k=1:size(mesh.ymin_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =0.0
              ny[i,kf] =-1.0

              ip = mesh.ymin_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              x = mesh.x[ip]
              y = mesh.y[ip]
              e = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,kf,1] = e
              exact_q[i,kf,2] = kx*e/c
              exact_q[i,kf,3] = ky*e/c
              flux_q[i,kf,1,1] = -((2*kx^2)*x+(2*kx*ky*(y-y0)-(2*kx^2)*x0-2*kx*c*t))*e
              flux_q[i,kf,2,1] = -((2*ky^2)*y+(2*kx*ky*(x-x0)-(2*ky^2)*y0-2*ky*c*t))*e
              flux_q[i,kf,1,2] = flux_q[i,kf,1,1]*kx/c
              flux_q[i,kf,2,2] = flux_q[i,kf,2,1]*kx/c
              flux_q[i,kf,1,3] = flux_q[i,kf,1,1]*ky/c
              flux_q[i,kf,2,3] = flux_q[i,kf,2,1]*ky/c
          end
      end
      disp = size(mesh.xmin_faces,2) + size(mesh.xmax_faces,2) + size(mesh.ymin_faces,2)
      for k=1:size(mesh.ymax_faces,2)
          for i=1:mesh.ngl
              kf = k+disp
              nx[i,kf] =0.0
              ny[i,kf] =1.0

              ip = mesh.ymax_faces[i,k]
              unl = nx[i,kf]*q[ip,2]+ny[i,kf]*q[ip,3]
              x = mesh.x[ip]
              y = mesh.y[ip]
              e = exp(- ((kx*(x - x0) + ky*(y - y0)-c*t)^2)/d2)
              unl = nx[i,k]*q[ip,2]+ny[i,k]*q[ip,3]
              exact_q[i,kf,1] = e
              exact_q[i,kf,2] = kx*e/c
              exact_q[i,kf,3] = ky*e/c
              flux_q[i,kf,1,1] = -((2*kx^2)*x+(2*kx*ky*(y-y0)-(2*kx^2)*x0-2*kx*c*t))*e
              flux_q[i,kf,2,1] = -((2*ky^2)*y+(2*kx*ky*(x-x0)-(2*ky^2)*y0-2*ky*c*t))*e
              flux_q[i,kf,1,2] = flux_q[i,kf,1,1]*kx/c
              flux_q[i,kf,2,2] = flux_q[i,kf,2,1]*kx/c
              flux_q[i,kf,1,3] = flux_q[i,kf,1,1]*ky/c
              flux_q[i,kf,2,3] = flux_q[i,kf,2,1]*ky/c
          end
      end
  #nx .= -nx
  #ny .= -ny
end


function build_custom_bcs!(t,mesh,q,qe,gradq,::NSD_3D,::DefaultBC,exact_q,flux_q,nx,ny,nz,nvars)
  npoin = mesh.npoin
  for ip = 1:mesh.npoin
      x = mesh.x[ip]
      y = mesh.y[ip]
      z = mesh.z[ip]
      nx[ip]=0.0
      ny[ip]=0.0
      nz[ip]=0.0
      if (AlmostEqual(x,mesh.xmin) && AlmostEqual(y,mesh.ymin) && AlmostEqual(z,mesh.zmin))
        nx[ip]=sqrt(1/3)
        ny[ip]=sqrt(1/3)
        nz[ip]=sqrt(1/3)
      elseif (AlmostEqual(x,mesh.xmin) && AlmostEqual(y,mesh.ymax),AlmostEqual(z,mesh.zmin))
        nx[ip] = sqrt(1/3)
        ny[ip] = -sqrt(1/3)
        nz[ip] = sqrt(1/3)
      elseif (AlmostEqual(x,mesh.xmax) && AlmostEqual(y,mesh.ymin),AlmostEqual(z,mesh.zmin))
        nx[ip] = -sqrt(1/3)
        ny[ip] = sqrt(1/3)
        nz[ip] = sqrt(1/3)
      elseif (AlmostEqual(x,mesh.xmax) && AlmostEqual(y,mesh.ymax),AlmostEqual(z,mesh.zmin))
        nx[ip] = -sqrt(1/3)
        ny[ip] = -sqrt(1/3)
        nz[ip] = sqrt(1/3)
      elseif (AlmostEqual(x,mesh.xmin) && AlmostEqual(y,mesh.ymin) && AlmostEqual(z,mesh.zmax))
        nx[ip]=sqrt(1/3)
        ny[ip]=sqrt(1/3) 
        nz[ip]=-sqrt(1/3)
      elseif (AlmostEqual(x,mesh.xmin) && AlmostEqual(y,mesh.ymax) && AlmostEqual(z,mesh.zmax))
        nx[ip] = sqrt(1/3)
        ny[ip] = -sqrt(1/3)
        nz[ip] = -sqrt(1/3)
      elseif (AlmostEqual(x,mesh.xmax) && AlmostEqual(y,mesh.ymin) && AlmostEqual(z,mesh.zmax))
        nx[ip] = -sqrt(1/3)
        ny[ip] = sqrt(1/3)
        nz[ip] = -sqrt(1/3)
      elseif (AlmostEqual(x,mesh.xmax) && AlmostEqual(y,mesh.ymax) && AlmostEqual(z,mesh.zmax))
        nx[ip] = -sqrt(1/3)
        ny[ip] = -sqrt(1/3)
        nz[ip] = sqrt(1/3)
      elseif (AlmostEqual(x,mesh.xmin) && AlmostEqual(y,mesh.ymin))
        nx[ip]=sqrt(0.5)
        ny[ip]=sqrt(0.5)
      elseif (AlmostEqual(x,mesh.xmin) && AlmostEqual(y,mesh.ymax))
        nx[ip] = sqrt(0.5)
        ny[ip] = -sqrt(0.5)
      elseif (AlmostEqual(x,mesh.xmax) && AlmostEqual(y,mesh.ymin))
        nx[ip] = -sqrt(0.5)
        ny[ip] = sqrt(0.5)
      elseif (AlmostEqual(x,mesh.xmax) && AlmostEqual(y,mesh.ymax))
        nx[ip] = -sqrt(0.5)
        ny[ip] = -sqrt(0.5)
      elseif (AlmostEqual(x,mesh.xmin) && AlmostEqual(z,mesh.zmin))
        nx[ip]=sqrt(0.5)
        nz[ip]=sqrt(0.5)
      elseif (AlmostEqual(x,mesh.xmin) && AlmostEqual(z,mesh.zmax))
        nx[ip] = sqrt(0.5)
        nz[ip] = -sqrt(0.5)
      elseif (AlmostEqual(x,mesh.xmax) && AlmostEqual(z,mesh.zmin))
        nx[ip] = -sqrt(0.5)
        nz[ip] = sqrt(0.5)
      elseif (AlmostEqual(x,mesh.xmax) && AlmostEqual(z,mesh.zmax))
        nx[ip] = -sqrt(0.5)
        nz[ip] = -sqrt(0.5)
      elseif (AlmostEqual(z,mesh.zmin) && AlmostEqual(y,mesh.ymin))
        nz[ip]=sqrt(0.5)
        ny[ip]=sqrt(0.5)
      elseif (AlmostEqual(z,mesh.zmin) && AlmostEqual(y,mesh.ymax))
        nz[ip] = sqrt(0.5)
        ny[ip] = -sqrt(0.5)
      elseif (AlmostEqual(z,mesh.zmax) && AlmostEqual(y,mesh.ymin))
        nz[ip] = -sqrt(0.5)
        ny[ip] = sqrt(0.5)
      elseif (AlmostEqual(z,mesh.zmax) && AlmostEqual(y,mesh.ymax))
        nz[ip] = -sqrt(0.5)
        ny[ip] = -sqrt(0.5) 
      elseif (AlmostEqual(x,mesh.xmin))
         nx[ip]=1.0
      elseif  (AlmostEqual(y,mesh.ymin))
         ny[ip] =1.0
      elseif  (AlmostEqual(x,mesh.xmax))
         nx[ip]=-1.0
      elseif  (AlmostEqual(y,mesh.ymax))
         ny[ip] =-1.0
      end
      if (nx[ip] != 0.0)
          qe[ip,2] = 0.0
      end
      if (ny[ip] != 0.0)
          qe[ip,3] = 0.0
      end
      if (nz[ip] != 0.0)
          qe[ip,4] = 0.0
      end
      unl = nx[ip]*q[ip,2]+ny[ip]*q[ip,3]+nz[ip]*q[ip,4]
      exact_q[ip,1] = -q[ip,1]
      exact_q[ip,2] = q[ip,2] - 2*unl*nx[ip]
      exact_q[ip,3] = q[ip,3] - 2*unl*ny[ip]
      exact_q[ip,4] = q[ip,4] - 2*unl*nz[ip]
      for var=1:nvars
        unl = nx[ip]*gradq[1,ip,var]+ny[ip]*gradq[2,ip,var]+nz*gradq[3,ip,var]
        flux_q[ip,1,var] = gradq[1,ip,var] - 2*gradq[1,ip,var]*unl
        flux_q[ip,2,var] = gradq[2,ip,var] - 2*gradq[2,ip,var]*unl
        flux_q[ip,3,var] = gradq[3,ip,var] - 2*gradq[3,ip,var]*unl
      end
  end

end
function build_custom_bcs(t,mesh,q,gradq,::NSD_2D,::DirichletExample)
  npoin = mesh.npoin
  exact_q = zeros(npoin,nvar)
  flux_q = zeros(npoin,2,nvar)
  for ip = 1:mesh.npoin
      x = mesh.x[ip]
      y = mesh.y[ip]
      exact_q[ip] = 0.0
      flux_q[ip,1,:] = gradq[ip,1,:]
      flux_q[ip,2,:] = gradq[ip,2,:]
      return exact_q,flux_q
  end
end

function build_custom_bcs(t,mesh,q,gradq,::NSD_3D,::DirichletExample)
  npoin = mesh.npoin
  exact_q = zeros(npoin,nvar)
  flux_q = zeros(npoin,3,nvar)
  for ip = 1:mesh.npoin
      x = mesh.x[ip]
      y = mesh.y[ip]
      exact_q[ip] = 0.0
      flux_q[ip,1,:] = gradq[ip,1,:]
      flux_q[ip,2,:] = gradq[ip,2,:]
      flux_q[ip,3,:] = gradq[ip,3,:]
      return exact_q,flux_q
  end
end
