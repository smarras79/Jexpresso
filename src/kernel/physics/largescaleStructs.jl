Base.@kwdef mutable struct St_LargeScaleTendencies{T <: AbstractFloat, dims1, backend}

    # WIP

    Rad_cool      = KernelAbstractions.zeros(backend,  T, dims1)
    T_adv         = KernelAbstractions.zeros(backend,  T, dims1)
    q_adv         = KernelAbstractions.zeros(backend,  T, dims1)

end

function allocate_LargeScaleTendencies(npoin, mesh, inputs, T, backend; lLST=false)

    if lLST
        dims1 = (Int64(npoin),1)
    else
        dims1 = (Int64(1))
    end

    LST = St_LargeScaleTendencies{T, dims1, backend}()
    if (lLST) 
        read_large_scale!(backend, inputs[:LST_files], LST, mesh)
    end

    return LST
end
