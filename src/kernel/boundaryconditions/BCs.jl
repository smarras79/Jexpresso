#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------


include("../../io/mod_inputs.jl")

include("../AbstractTypes.jl")
include("../basis/basis_structs.jl")
include("../infrastructure/element_matrices.jl")
include("../infrastructure/Kopriva_functions.jl")
include("../infrastructure/2D_3D_structures.jl")
include("../mesh/metric_terms.jl")
include("../mesh/mesh.jl")



function apply_boundary_conditions!(qp,mesh,inputs,::NSD_2D)
   #Periodic boundaries
   if (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="periodic" || haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="periodic")
      inputs[:xmax_bc] = "periodic"
      inputs[:xmin_bc] = "periodic"
      # B.C.: periodic
      for I in keys(mesh.xperiodicity)
         ip = mesh.xperiodicity[I]
         qp.qn[I,:] .= 0.5*(qp.qn[I,:] .+ qp.qn[ip,:])
         qp.qn[ip,:] .= qp.qn[I,:]
      end
    end 
    if (haskey(inputs, :ymin_bc) && inputs[:ymin_bc]=="periodic") || (haskey(inputs, :ymax_bc) && inputs[:ymax_bc]=="periodic")
      inputs[:ymax_bc] = "periodic"
      inputs[:ymin_bc] = "periodic"
      for I in keys(mesh.yperiodicity)
         ip = mesh.yperiodicity[I]
         qp.qn[I,:] .= 0.5*(qp.qn[I,:] .+ qp.qn[ip,:])
         qp.qn[ip,:] .= qp.qn[I,:]
      end
    end

   #Dirichlet boundaries
   if (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="dirichlet") 
      #Left boundary
      exact = inputs[:bc_exact_xmin] 
      for I=1:size(mesh.bc_xmin,1)
         ip = mesh.bc_xmin[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
   if (haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="dirichlet") 
            #Right boundary
      exact = inputs[:bc_exact_xmax]
      for I=1:size(mesh.bc_xmax,1)
         ip = mesh.bc_xmax[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
   if (haskey(inputs, :ymin_bc) && inputs[:ymin_bc]=="dirichlet")
         #bottom boundary
      exact = inputs[:bc_exact_ymin]
      for I=1:size(mesh.bc_ymin,1)
         ip = mesh.bc_ymin[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
   if (haskey(inputs, :ymax_bc) && inputs[:ymax_bc]=="dirichlet") 
            #top boundary
      exact = inputs[:bc_exact_ymax]
      for I=1:size(mesh.bc_ymax,1)
         ip = mesh.bc_ymax[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
end

function apply_boundary_conditions!(qp,mesh,inputs,::NSD_3D)

    if (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="periodic" || haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="periodic")
      inputs[:xmax_bc] = "periodic"
      inputs[:xmin_bc] = "periodic"
      # B.C.: periodic
      for I in keys(mesh.xperiodicity)
         ip = mesh.xperiodicity[I]
         qp.qn[I,:] .= 0.5*(qp.qn[I,:] .+ qp.qn[ip,:])
         qp.qn[ip,:] .= qp.qn[I,:]
      end
    end
    if (haskey(inputs, :ymin_bc) && inputs[:ymin_bc]=="periodic") || (haskey(inputs, :ymax_bc) && inputs[:ymax_bc]=="periodic")
      inputs[:ymax_bc] = "periodic"
      inputs[:ymin_bc] = "periodic"
      for I in keys(mesh.yperiodicity)
         ip = mesh.yperiodicity[I]
         qp.qn[I,:] .= 0.5*(qp.qn[I,:] .+ qp.qn[ip,:])
         qp.qn[ip,:] .= qp.qn[I,:]
      end
    end
    if (haskey(inputs, :zmin_bc) && inputs[:zmin_bc]=="periodic") || (haskey(inputs, :zmax_bc) && inputs[:zmax_bc]=="periodic")
      inputs[:zmax_bc] = "periodic"
      inputs[:zmin_bc] = "periodic"
      for I in keys(mesh.zperiodicity)
         ip = mesh.zperiodicity[I]
         qp.qn[I,:] .= 0.5*(qp.qn[I,:] .+ qp.qn[ip,:])
         qp.qn[ip,:] .= qp.qn[I,:]
      end
    end

   #Dirichlet boundaries
   if (haskey(inputs, :xmin_bc) && inputs[:xmin_bc]=="dirichlet")
      #Left boundary
      exact = inputs[:bc_exact_xmin]
      for I=1:size(mesh.bc_xmin,1)
         ip = mesh.bc_xmin[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
   if (haskey(inputs, :xmax_bc) && inputs[:xmax_bc]=="dirichlet") 
            #Right boundary
      exact = inputs[:bc_exact_xmax]
      for I=1:size(mesh.bc_xmax,1)
         ip = mesh.bc_xmax[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
   if (haskey(inputs, :ymin_bc) && inputs[:ymin_bc]=="dirichlet")
         #bottom boundary
      exact = inputs[:bc_exact_ymin]
      for I=1:size(mesh.bc_ymin,1)
         ip = mesh.bc_ymin[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
   if (haskey(inputs, :ymax_bc) && inputs[:ymax_bc]=="dirichlet")
            #top boundary
      exact = inputs[:bc_exact_ymax]
      for I=1:size(mesh.bc_ymax,1)
         ip = mesh.bc_ymax[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
   if (haskey(inputs, :zmin_bc) && inputs[:zmin_bc]=="dirichlet")
         #bottom boundary
      exact = inputs[:bc_exact_zmin]
      for I=1:size(mesh.bc_zmin,1)
         ip = mesh.bc_zmin[I]
         qp.qn[ip,:] .= exact[:]
      end
   end
   if (haskey(inputs, :zmax_bc) && inputs[:zmax_bc]=="dirichlet")
            #top boundary
      exact = inputs[:bc_exact_zmax]
      for I=1:size(mesh.bc_zmax,1)
         ip = mesh.bc_zmax[I]
         qp.qn[ip,:] .= exact[:]
      end
   end

end
