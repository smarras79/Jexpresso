function convert_mesh_arrays!(mesh, backend)

    aux = KernelAbstractions.allocate(backend, TFloat, mesh.npoin)
    KernelAbstractions.copyto!(backend, aux, mesh.x)
    mesh.x = KernelAbstractions.allocate(backend, TFloat, mesh.npoin)
    mesh.x .= aux

    aux = KernelAbstractions.allocate(backend, TFloat, mesh.npoin)
    KernelAbstractions.copyto!(backend, aux, mesh.y)
    mesh.y = KernelAbstractions.allocate(backend, TFloat, mesh.npoin)
    mesh.y .= aux

    aux = KernelAbstractions.allocate(backend, TFloat, mesh.npoin)
    KernelAbstractions.copyto!(backend, aux, mesh.z)
    mesh.z = KernelAbstractions.allocate(backend, TFloat, mesh.npoin)
    mesh.z .= aux

    aux = KernelAbstractions.allocate(backend, TInt, mesh.nelem, mesh.ngl, mesh.ngl)
    KernelAbstractions.copyto!(backend, aux, mesh.connijk)
    mesh.connijk = KernelAbstractions.allocate(backend, TInt, mesh.nelem, mesh.ngl, mesh.ngl)
    mesh.connijk .= aux

end

function convert_mesh_arrays_to_cpu!(mesh)

    aux = KernelAbstractions.allocate(CPU(), TFloat, mesh.npoin)
    KernelAbstractions.copyto!(CPU(), aux, mesh.x)
    mesh.x = KernelAbstractions.allocate(CPU(), TFloat, mesh.npoin)
    mesh.x .= aux

    aux = KernelAbstractions.allocate(CPU(), TFloat, mesh.npoin)
    KernelAbstractions.copyto!(CPU(), aux, mesh.y)
    mesh.y = KernelAbstractions.allocate(CPU(), TFloat, mesh.npoin)
    mesh.y .= aux

    aux = KernelAbstractions.allocate(CPU(), TFloat, mesh.npoin)
    KernelAbstractions.copyto!(CPU(), aux, mesh.z)
    mesh.z = KernelAbstractions.allocate(CPU(), TFloat, mesh.npoin)
    mesh.z .= aux

    aux = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem, mesh.ngl, mesh.ngl)
    KernelAbstractions.copyto!(CPU(), aux, mesh.connijk)
    mesh.connijk = KernelAbstractions.allocate(CPU(), TInt, mesh.nelem, mesh.ngl, mesh.ngl)
    mesh.connijk .= aux
end

