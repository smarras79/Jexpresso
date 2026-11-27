Base.@kwdef mutable struct St_Wall_model{T <: AbstractFloat, dims1, dims2, dims3, backend, VT}

    τ_f::VT = KernelAbstractions.zeros(backend,  T, dims1)
    wθ::VT  = KernelAbstractions.zeros(backend,  T, dims2)
    wqv::VT = KernelAbstractions.zeros(backend,  T, dims3)

end

# Whether any boundary face of the (whole, MPI-partitioned) 3D mesh is tagged "MOST", the tag
# the Monin-Obukhov wall model acts on. The wall-model arrays below, and the LES statistics of the
# MOST surface, exist only then. Agreed on by all ranks (both uses involve collectives or must be
# consistent across ranks). 1D/2D: false (the 2D MOST branch uses local work arrays).
function mesh_has_most_tag(mesh)
    local_has = mesh.SD == NSD_3D() && any(==("MOST"), mesh.bdy_face_type)
    return MPI.Allreduce(local_has, MPI.LOR, get_mpi_comm())
end

function allocate_Wall_model(nface, ngl, T, backend; lwall_model=false, lmoist=false)
    # Fixed: Only allocate arrays if wall model is enabled to save memory
    if (lwall_model == true) && (lmoist == true)
        dims1 = (nface, ngl, ngl, 3)
        dims2 = (nface, ngl, ngl, 1)
        dims3 = (nface, ngl, ngl, 1)
    elseif (lwall_model == true) && (lmoist == false)
        dims1 = (nface, ngl, ngl, 3)
        dims2 = (nface, ngl, ngl, 1)
        dims3 = (0, 0, 0, 0)
    else
        # Allocate minimal empty arrays when wall model is disabled
        dims1 = (0, 0, 0, 0)
        dims2 = (0, 0, 0, 0)
        dims3 = (0, 0, 0, 0)
    end

    VT = typeof(KernelAbstractions.zeros(backend, T, dims1))
    wm = St_Wall_model{T, dims1, dims2, dims3, backend, VT}()
    return wm
end
