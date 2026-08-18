#
# mesh.x/y/z are VIEWS of rows of mesh.coords (see meshStructs.jl), so staging
# coords stages them. The three per-array copy blocks that used to be here moved
# the same numbers a second time and are gone; `mesh.x = …` would now raise.
#
function convert_mesh_arrays!(::NSD_1D, mesh, backend, inputs)
    mesh.coords = convert_to_typed_array(mesh.coords, TFloat)

    #
    # coords (NSD_MAX x npoin) is THE coordinate storage — mesh.x/y/z are views of
    # its rows. The device kernels take a node's coordinates as
    # @view(coords[:, ip]), one cache line per node, which is the point of the
    # layout.
    #
    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.coords))
    KernelAbstractions.copyto!(backend, aux, mesh.coords)
    mesh.coords = KernelAbstractions.allocate(backend, TFloat, size(mesh.coords))
    mesh.coords .= aux

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
    mesh.coords = convert_to_typed_array(mesh.coords, TFloat)

    #
    # coords (NSD_MAX x npoin) is THE coordinate storage — mesh.x/y/z are views of
    # its rows. The device kernels take a node's coordinates as
    # @view(coords[:, ip]), one cache line per node, which is the point of the
    # layout.
    #
    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.coords))
    KernelAbstractions.copyto!(backend, aux, mesh.coords)
    mesh.coords = KernelAbstractions.allocate(backend, TFloat, size(mesh.coords))
    mesh.coords .= aux

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
    mesh.coords = convert_to_typed_array(mesh.coords, TFloat)
    npoin = size(getfield(mesh, :coords), 2)
    #
    # coords (NSD_MAX x npoin) is THE coordinate storage — mesh.x/y/z are views of
    # its rows. The device kernels take a node's coordinates as
    # @view(coords[:, ip]), one cache line per node, which is the point of the
    # layout.
    #
    aux = KernelAbstractions.allocate(backend, TFloat, size(mesh.coords))
    KernelAbstractions.copyto!(backend, aux, mesh.coords)
    mesh.coords = KernelAbstractions.allocate(backend, TFloat, size(mesh.coords))
    mesh.coords .= aux

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

    # coords comes home with x/y/z — see convert_mesh_arrays!.
    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.coords))
    KernelAbstractions.copyto!(CPU(), aux, mesh.coords)
    mesh.coords = KernelAbstractions.allocate(CPU(), Float64, size(mesh.coords))
    mesh.coords .= aux

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

    # coords comes home with x/y/z — see convert_mesh_arrays!.
    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.coords))
    KernelAbstractions.copyto!(CPU(), aux, mesh.coords)
    mesh.coords = KernelAbstractions.allocate(CPU(), Float64, size(mesh.coords))
    mesh.coords .= aux

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

    npoin = size(getfield(mesh, :coords), 2)
    # coords comes home with x/y/z — see convert_mesh_arrays!.
    aux = KernelAbstractions.allocate(CPU(), Float64, size(mesh.coords))
    KernelAbstractions.copyto!(CPU(), aux, mesh.coords)
    mesh.coords = KernelAbstractions.allocate(CPU(), Float64, size(mesh.coords))
    mesh.coords .= aux

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
