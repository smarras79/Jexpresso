using LinearSolve
using SnoopCompile
using WriteVTK
import SciMLBase

include("./plotting/jeplots.jl")

abstract type AbstractOutFormat end
struct PNG <: AbstractOutFormat end
struct ASCII <: AbstractOutFormat end 

#----------------------------------------------------------------------------------------------------------------------------------------------
# ∂q/∂t = RHS -> q(x,t)
#----------------------------------------------------------------------------------------------------------------------------------------------
# PNG
function write_output(sol::ODESolution, SD::NSD_3D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1) nothing end
function write_output(sol::ODESolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    
    if inputs[:lplot_surf3d]
        for iout = 1:inputs[:ndiagnostics_outputs]
            title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
            plot_surf3d(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR; iout=iout, nvar=nvar, smoothing_factor=inputs[:smoothing_factor])
        end
    else
        for iout = 1:inputs[:ndiagnostics_outputs]
            title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
            plot_triangulation(SD, mesh, sol.u[iout][:], title,  OUTPUT_DIR; iout=iout, nvar=nvar)
        end
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE"))
end
function write_output(sol::ODESolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)
    
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  "))
    for iout = 1:inputs[:ndiagnostics_outputs]
        title = string("sol.u at time ", sol.t[iout])
        plot_results(SD, mesh, sol.u[iout][:], title, OUTPUT_DIR; iout=iout, nvar=nvar)
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE ") )
end

# ASCII
function write_output(sol::ODESolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1)
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    for iout = 1:inputs[:ndiagnostics_outputs]
        title = string("sol.u at time ", sol.t[iout])

        for ivar=1:nvar
            fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".dat")
            open(fout_name, "w") do f
                idx = (ivar - 1)*mesh.npoin
                for ip = (ivar-1)*mesh.npoin+1:ivar*mesh.npoin
                    @printf(f, " %f \n", sol.u[iout][ip])
                end
            end
        end
        #write_ascii(SD, mesh.x, sol.u[iout][:], title, OUTPUT_DIR; iout=iout, nvar=nvar)
    end
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  DONE ") )
end
function write_output(sol::ODESolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1)
    
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    for iout = 1: inputs[:ndiagnostics_outputs]
        #Write out data at final timestep
	fname = @sprintf "it-%d.dat" iout
    	open(string(OUTPUT_DIR, "/", fname), "w") do f
            for ip = 1:mesh.npoin
                @printf(f, " %d %.6f %.6f %.6f \n", ip, mesh.x[ip], mesh.y[ip], sol.u[iout][ip])
            end
        end #f
    end    
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  DONE ") ) 
end




#----------------------------------------------------------------------------------------------------------------------------------------------
# Aq = b -> q(x)
#----------------------------------------------------------------------------------------------------------------------------------------------
# PNG
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::PNG; nvar=1)

    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  ") )
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    
    if inputs[:lplot_surf3d]
        plot_surf3d(SD, mesh, sol.u, title, OUTPUT_DIR; iout=1, nvar=1, smoothing_factor=inputs[:smoothing_factor])
    else
        plot_triangulation(SD, mesh, sol.u, title, OUTPUT_DIR;)
    end
    println(string(" # Writing output to PNG file:", OUTPUT_DIR, "*.png ...  DONE") )
end

# ASCII
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_1D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1) nothing end
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_3D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1) nothing end
function write_output(sol::SciMLBase.LinearSolution, SD::NSD_2D, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict, outformat::ASCII; nvar=1)   
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  ") )
    title = @sprintf "Solution to ∇⋅∇(q) = f"
    fname = @sprintf "Axb.dat"
    open(string(OUTPUT_DIR, "/", fname), "w") do f
        for ip = 1:length(sol.u)
            @printf(f, " %d %.6f %.6f %.6f \n", ip, mesh.x[ip], mesh.y[ip], sol.u[ip])
        end #f
    end
    println(string(" # Writing output to ASCII file:", OUTPUT_DIR, "*.dat ...  DONE") )
end

#VTK:
function write_vtk(sol::Array, SD, mesh::St_mesh, OUTPUT_DIR::String, inputs::Dict)
    nothing
    #println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ...  ") )
    #cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:mesh.npoin]
    #vtk_grid(string(OUTPUT_DIR, "qsolution"), mesh.x, mesh.y, mesh.z, cells) do vtk
    #   vtk["q", VTKPointData()] = sol[:,1]
    #end
    #println(string(" # Writing output to VTK file:", OUTPUT_DIR, "*.vtu ...  DONE") )
end

#---------------------------------------------------------------------!
# This function creates cells array for VTK using mesh.conn
# 
# From numa Michal A. Kopera on 01/2015
#---------------------------------------------------------------------!
function create_cells(mesh::St_mesh, cells,ncells, mesh::St_mesh)

    
    cells = zeros(Int64, ncells)
    ic=0

    if (mesh.SD === NSD_2D())
        for e=1,nelem
            #construct cells in each element
            for k=1;maximum(nglz-1,1)
                for j=1:maximum(ngly-1,1)
                    for i=1:maximum(nglx-1,1)
                        ii=min(i+1,nglx)
                        jj=min(j+1,ngly)
                        kk=min(k+1,nglz)
                        ic=ic+1

                        if(nglx == 1) 
                            cells(1,ic)= intma( i, j, k,e) -1
                            cells(2,ic)= intma( i,jj, k,e) -1
                            cells(3,ic)= intma( i,jj,kk,e) -1
                            cells(4,ic)= intma( i, j,kk,e) -1
                        elseif(ngly == 1) 
                            cells(1,ic)= intma( i, j, k,e) -1
                            cells(2,ic)= intma(ii, j, k,e) -1
                            cells(3,ic)= intma(ii, j,kk,e) -1
                            cells(4,ic)= intma( i, j,kk,e) -1
                        elseif(nglz == 1) 
                            cells(1,ic)= intma( i, j, k,e) -1
                            cells(2,ic)= intma(ii, j, k,e) -1
                            cells(3,ic)= intma(ii,jj, k,e) -1
                            cells(4,ic)= intma( i,jj, k,e) -1
                        else
                            cells(1,ic)= intma( i, j, k,e) -1
                            cells(2,ic)= intma(ii, j, k,e) -1
                            cells(3,ic)= intma(ii,jj, k,e) -1
                            cells(4,ic)= intma( i,jj, k,e) -1
                            cells(5,ic)= intma( i, j,kk,e) -1
                            cells(6,ic)= intma(ii, j,kk,e) -1
                            cells(7,ic)= intma(ii,jj,kk,e) -1
                            cells(8,ic)= intma( i,jj,kk,e) -1
                            endif

                        end !i
                    end !j
                end 
            end  !ie
        end
    end
    
end

function outvtk(mesh::St_mesh, q, q_ref, fname; nvar)
  
    #local arrays
  real u, v, pp
  integer ivar, ifactor, i, j, e, ie, ncells, nsize
  integer i1, i2, i3, i4
  real rfactor
  
  #Open VTK file
  open(1,file=fname)
  
  write(1,'(a)')"# vtk DataFile Version 2.0"
  write(1,'(a)')'Numa2dCG Data'
  write(1,'(a)')'ASCII'
  write(1,'(a)')'DATASET UNSTRUCTURED_GRID'
  
  #Write out point field
  write(1,*)'POINTS',npoin,'float'
  for i=1,npoin
     #write(1,'(3(e12.4,1x))')coord(1,i), coord(2,i), 0.0
     write(1,*) coord(1,i), coord(2,i), 0.0
  end
  
  #Write out cell field
  ncells=nelem*(ngl-1)*(ngl-1)
  nsize=ncells*5
  write(1,*)"CELLS",ncells,nsize
  for ie =1,nelem
     for i = 1,ngl-1
        for j = 1,ngl-1
            write(1,*)'4', &
                (mesh.connijk(i,j,ie)-1), &
                (mesh.connijk(i+1,j,ie)-1),  &
                (mesh.connijk(i+1,j+1,ie)-1),  &
                (mesh.connijk(i,j+1,ie)-1)
        end
     end
  end
 
  write(1,*)'CELL_TYPES',ncells
  for i=1,ncells
     write(1,*)'9'
  end
  
  write(1,*)'CELL_DATA', ncells
  write(1,*)'POINT_DATA', npoin
  
  write(1,*)'SCALARS dRho foruble'
  write(1,*)'LOOKUP_TABLE default'
  for i=1,npoin
     #write(1,'((1(e12.2,1x)))')q(1,i)
     write(1,*) q(1,i)
  end

  write(1,*)'SCALARS Rho foruble'
  write(1,*)'LOOKUP_TABLE default'
  for i=1,npoin
     #write(1,'((1(e12.2,1x)))')q(1,i)
     write(1,*) q(1,i) + q_ref(1,i)
  end
  
  write(1,*)'SCALARS Uvelo foruble'
  write(1,*)'LOOKUP_TABLE default'
  for i=1,npoin
     #write(1,'((1(e12.2,1x)))')q(2,i)
     write(1,*) q(2,i)
  end

  write(1,*)'SCALARS Vvelo foruble'
  write(1,*)'LOOKUP_TABLE default'
  for i=1,npoin
     #write(1,'((1(e12.2,1x)))')q(3,i)
     write(1,*) q(3,i)
  end

  write(1,*)'SCALARS dTheta foruble'
  write(1,*)'LOOKUP_TABLE default'
  for i=1,npoin
     #write(1,'((1(e12.2,1x)))')q(4,i)
     write(1,*) q(4,i)
  end

  write(1,*)'SCALARS Theta foruble'
  write(1,*)'LOOKUP_TABLE default'
  for i=1,npoin
     #write(1,'((1(e12.2,1x)))')q(4,i)
     write(1,*) q(4,i) + q_ref(4,i)
  end
  
 
  write(1,*)'VECTORS Velocity foruble'
  for i=1,npoin
     #write(1,'((3(e12.4,1x)),a)')q(2,i),q(3,i),0
     write(1,*) q(2,i), q(3,i), 0.0
  end
    
end subroutine outvtk_time

