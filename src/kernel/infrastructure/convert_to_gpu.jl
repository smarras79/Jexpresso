function convert_mesh_arrays!(::NSD_1D, mesh, backend, inputs)
    # Ensure mesh.x, mesh.y, mesh.z are of type TFloat (e.g., Float32)
    mesh.x = convert_to_typed_array(mesh.x, TFloat)
    mesh.y = convert_to_typed_array(mesh.y, TFloat)
    mesh.z = convert_to_typed_array(mesh.z, TFloat)

    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.x))
    KernelAbstractions.copyto!(backend, aux, mesh.x)
    mesh.x = KernelAbstractions.allocate(backend, TFloat, size(mesh.x))
    mesh.x .= aux

    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.y))
    KernelAbstractions.copyto!(backend, aux, mesh.y)
    mesh.y = KernelAbstractions.allocate(backend, TFloat, size(mesh.y))
    mesh.y .= aux

    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.z))
    KernelAbstractions.copyto!(backend, aux, mesh.z)
    mesh.z = KernelAbstractions.allocate(backend, TFloat, size(mesh.z))
    mesh.z .= aux

    aux = KernelAbstractions.allocate(backend, TInt, mesh.nelem, mesh.ngl, 1)
    KernelAbstractions.copyto!(backend, aux, mesh.connijk)
    mesh.connijk = KernelAbstractions.allocate(backend, TInt, mesh.nelem, mesh.ngl, 1)
    mesh.connijk .= aux

    if (inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])
        aux = KernelAbstractions.allocate(backend, TInt, mesh.nelem_semi_inf, mesh.ngr, 1)
        KernelAbstractions.copyto!(backend, aux, mesh.connijk_lag)
        mesh.connijk_lag = KernelAbstractions.allocate(backend, TInt, mesh.nelem_semi_inf, mesh.ngr, 1)
        mesh.connijk_lag .= aux
    end

end

function convert_mesh_arrays!(::NSD_2D, mesh, backend, inputs)
    # Ensure mesh.x, mesh.y, mesh.z are of type TFloat (e.g., Float32)
    mesh.x = convert_to_typed_array(mesh.x, TFloat)
    mesh.y = convert_to_typed_array(mesh.y, TFloat)
    mesh.z = convert_to_typed_array(mesh.z, TFloat)

    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.x))
    KernelAbstractions.copyto!(backend, aux, mesh.x)
    mesh.x = KernelAbstractions.allocate(backend, TFloat, size(mesh.x))
    mesh.x .= aux

    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.y))
    KernelAbstractions.copyto!(backend, aux, mesh.y)
    mesh.y = KernelAbstractions.allocate(backend, TFloat, size(mesh.y))
    mesh.y .= aux

    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.z))
    KernelAbstractions.copyto!(backend, aux, mesh.z)
    mesh.z = KernelAbstractions.allocate(backend, TFloat, size(mesh.z))
    mesh.z .= aux

    aux = KernelAbstractions.allocate(backend, TInt, mesh.nelem, mesh.ngl, mesh.ngl)
    KernelAbstractions.copyto!(backend, aux, mesh.connijk)
    mesh.connijk = KernelAbstractions.allocate(backend, TInt, mesh.nelem, mesh.ngl, mesh.ngl)
    mesh.connijk .= aux

    if ("Laguerre" in mesh.bdy_edge_type)
        aux = KernelAbstractions.allocate(backend, TInt, mesh.nelem_semi_inf, mesh.ngl, mesh.ngr)
        KernelAbstractions.copyto!(backend, aux, mesh.connijk_lag)
        mesh.connijk_lag = KernelAbstractions.allocate(backend, TInt, mesh.nelem_semi_inf, mesh.ngl, mesh.ngr)
        mesh.connijk_lag .= aux
    end

    aux = KernelAbstractions.allocate(backend, TInt, mesh.nedges_bdy, mesh.ngl)
    KernelAbstractions.copyto!(backend, aux, mesh.poin_in_bdy_edge)
    mesh.poin_in_bdy_edge = KernelAbstractions.allocate(backend, TInt, mesh.nedges_bdy, mesh.ngl)
    mesh.poin_in_bdy_edge .= aux
end

function convert_mesh_arrays!(::NSD_3D, mesh, backend, inputs)
    # Ensure mesh.x, mesh.y, mesh.z are of type TFloat (e.g., Float32)
    mesh.x = convert_to_typed_array(mesh.x, TFloat)
    mesh.y = convert_to_typed_array(mesh.y, TFloat)
    mesh.z = convert_to_typed_array(mesh.z, TFloat)
    npoin = size(mesh.x, 1)
    aux = KernelAbstractions.allocate(backend, TFloat, npoin)
    KernelAbstractions.copyto!(backend, aux, mesh.x)
    mesh.x = KernelAbstractions.allocate(backend, TFloat, npoin)
    mesh.x .= aux

    aux = KernelAbstractions.allocate(backend, TFloat, npoin)
    KernelAbstractions.copyto!(backend, aux, mesh.y)
    mesh.y = KernelAbstractions.allocate(backend, TFloat, npoin)
    mesh.y .= aux

    aux = KernelAbstractions.allocate(backend, TFloat, npoin)
    KernelAbstractions.copyto!(backend, aux, mesh.z)
    mesh.z = KernelAbstractions.allocate(backend, TFloat, npoin)
    mesh.z .= aux

    aux = KernelAbstractions.allocate(backend, TInt, mesh.nelem, mesh.ngl, mesh.ngl, mesh.ngl)
    KernelAbstractions.copyto!(backend, aux, mesh.connijk)
    mesh.connijk = KernelAbstractions.allocate(backend, TInt, mesh.nelem, mesh.ngl, mesh.ngl, mesh.ngl)
    mesh.connijk .= aux

    aux = KernelAbstractions.allocate(backend, TInt, mesh.nfaces_bdy, mesh.ngl, mesh.ngl)
    KernelAbstractions.copyto!(backend, aux, mesh.poin_in_bdy_face)
    mesh.poin_in_bdy_face = KernelAbstractions.allocate(backend, TInt, mesh.nfaces_bdy, mesh.ngl, mesh.ngl)
    mesh.poin_in_bdy_face .= aux
end

function convert_mesh_arrays_to_cpu!(::NSD_1D, mesh, inputs)

    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.x))
    KernelAbstractions.copyto!(CPU(), aux, mesh.x)
    mesh.x = KernelAbstractions.allocate(CPU(), Float64, size(mesh.x))
    mesh.x .= aux

    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.y))
    KernelAbstractions.copyto!(CPU(), aux, mesh.y)
    mesh.y = KernelAbstractions.allocate(CPU(), Float64, size(mesh.y))
    mesh.y .= aux

    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.z))
    KernelAbstractions.copyto!(CPU(), aux, mesh.z)
    mesh.z = KernelAbstractions.allocate(CPU(), Float64, size(mesh.z))
    mesh.z .= aux

    aux = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem, mesh.ngl, 1)
    KernelAbstractions.copyto!(CPU(), aux, mesh.connijk)
    mesh.connijk = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem, mesh.ngl, 1)
    mesh.connijk .= aux

    if (inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])
        aux = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem_semi_inf, mesh.ngr, 1)
        KernelAbstractions.copyto!(CPU(), aux, mesh.connijk_lag)
        mesh.connijk_lag = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem_semi_inf, mesh.ngr, 1)
        mesh.connijk_lag .= aux
    end

end

function convert_mesh_arrays_to_cpu!(::NSD_2D, mesh, inputs)

    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.x))
    KernelAbstractions.copyto!(CPU(), aux, mesh.x)
    mesh.x = KernelAbstractions.allocate(CPU(), Float64, size(mesh.x))
    mesh.x .= aux

    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.y))
    KernelAbstractions.copyto!(CPU(), aux, mesh.y)
    mesh.y = KernelAbstractions.allocate(CPU(), Float64, size(mesh.y))
    mesh.y .= aux

    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.z))
    KernelAbstractions.copyto!(CPU(), aux, mesh.z)
    mesh.z = KernelAbstractions.allocate(CPU(), Float64, size(mesh.z))
    mesh.z .= aux

    aux = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem, mesh.ngl, mesh.ngl)
    KernelAbstractions.copyto!(CPU(), aux, mesh.connijk)
    mesh.connijk = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem, mesh.ngl, mesh.ngl)
    mesh.connijk .= aux

    if ("Laguerre" in mesh.bdy_edge_type)
        aux = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem_semi_inf, mesh.ngl, mesh.ngr)
        KernelAbstractions.copyto!(CPU(), aux, mesh.connijk_lag)
        mesh.connijk_lag = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem_semi_inf, mesh.ngl, mesh.ngr)
        mesh.connijk_lag .= aux
    end

    aux = KernelAbstractions.allocate(CPU(), TInt, mesh.nedges_bdy, mesh.ngl)
    KernelAbstractions.copyto!(CPU(), aux, mesh.poin_in_bdy_edge)
    mesh.poin_in_bdy_edge = KernelAbstractions.allocate(CPU(), TInt, mesh.nedges_bdy, mesh.ngl)
    mesh.poin_in_bdy_edge .= aux
end

function convert_mesh_arrays_to_cpu!(::NSD_3D, mesh, inputs)

    npoin = size(mesh.x,1)
    aux = KernelAbstractions.allocate(CPU(), TFloat, npoin)
    KernelAbstractions.copyto!(CPU(), aux, mesh.x)
    mesh.x = KernelAbstractions.allocate(CPU(), TFloat, npoin)
    mesh.x .= aux

    aux = KernelAbstractions.allocate(CPU(), TFloat, npoin)
    KernelAbstractions.copyto!(CPU(), aux, mesh.y)
    mesh.y = KernelAbstractions.allocate(CPU(), TFloat, npoin)
    mesh.y .= aux

    aux = KernelAbstractions.allocate(CPU(), TFloat, npoin)
    KernelAbstractions.copyto!(CPU(), aux, mesh.z)
    mesh.z = KernelAbstractions.allocate(CPU(), TFloat, npoin)
    mesh.z .= aux

    aux = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem, mesh.ngl, mesh.ngl, mesh.ngl)
    KernelAbstractions.copyto!(CPU(), aux, mesh.connijk)
    mesh.connijk = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem, mesh.ngl, mesh.ngl, mesh.ngl)
    mesh.connijk .= aux

    aux = KernelAbstractions.allocate(CPU(), TInt, mesh.nfaces_bdy, mesh.ngl, mesh.ngl)
    KernelAbstractions.copyto!(CPU(), aux, mesh.poin_in_bdy_face)
    mesh.poin_in_bdy_face = KernelAbstractions.allocate(CPU(), TInt, mesh.nfaces_bdy, mesh.ngl, mesh.ngl)
    mesh.poin_in_bdy_face .= aux
end

# Helper function to convert arrays to the desired type
function convert_to_typed_array(arr, T)
    return T.(arr)
end
