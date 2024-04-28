
abstract type AbstractSpaceDimensions end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

Base.@kwdef mutable struct St_a{TFloat <: AbstractFloat}
    x::Array{TFloat} = zeros(TFloat, 0)
end

function allocate_operator_arrays(SD::NSD_1D, nelem, ngl, TFloat; neqs=1)
    
    #a = zeros(TFloat, nelem, ngl, neqs)
    a = St_a{TFloat}(x = zeros(TFloat, nelem, ngl))

    @info""
    @info " 1D"
    @info size(a.x)
    return a
end

function allocate_operator_arrays(SD::NSD_2D, nelem, ngl, TFloat; neqs=1)
    
    #a = zeros(TFloat, nelem, ngl, ngl, neqs)
    a = St_a{TFloat}(x = zeros(TFloat, nelem, ngl, ngl))
    @info ""
    @info " 2D"
    @info size(a.x)

    return a
end


function allocate_operator_arrays(SD::NSD_3D, nelem, ngl, TFloat; neqs=1)
    
    #a = zeros(TFloat, nelem, ngl, ngl, neqs)
    a = St_a{TFloat}(x = zeros(TFloat, nelem, ngl, ngl, ngl))
    @info ""
    @info " 3D"
    @info size(a.x)

    return a
end


function myfun!(params)
    for i=1:params.nelem
        for j=1:params.ngl
            for k=1:params.ngl
                for l=1:params.ngl
                    params.a.x[i,j,k,l] = 25.0
                end
            end
        end
    end
end

function myfun!(x, nelem, ngl)
    for i=1:nelem
        for j=1:ngl
            for k=1:ngl
                for l=1:ngl
                    x[i,j,k,l] = 25.0
                end
            end
        end
    end
end

function main()

    nelem = 10
    ngl   =  5
    a = allocate_operator_arrays(NSD_1D(), nelem, ngl, Float64; neqs=1)

    @info ""
    @info " main 1D"
    @info size(a.x)
    
    a = allocate_operator_arrays(NSD_2D(), nelem, ngl, Float64; neqs=1)

    @info ""
    @info " main 2D"
    @info size(a.x)
     
    a = allocate_operator_arrays(NSD_3D(), nelem, ngl, Float64; neqs=1)

    @info ""
    @info " main 3D"
    @info size(a.x)

    params = (a=a, nelem=nelem, ngl=ngl)
    @info params.nelem
    
    @info " 111"
    @time myfun!(a.x, nelem, ngl)

    @info " 222"
    @time myfun!(params)
    
    
end

main()
