
abstract type AbstractSpaceDimensions end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

Base.@kwdef mutable struct St_a{TFloat <: AbstractFloat}
    x::Array{TFloat} = zeros(TFloat, 0)
end


Base.@kwdef mutable struct St_flux{TFloat <: AbstractFloat}
    F::Array{TFloat} = zeros(TFloat, 0)
    G::Array{TFloat} = zeros(TFloat, 0)
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
    
    a    = St_a{TFloat}(x = zeros(TFloat, nelem, ngl, ngl, ngl))
    flux = St_flux{TFloat}(F = zeros(TFloat, nelem, ngl, ngl, ngl),
                        G = zeros(TFloat, nelem, ngl, ngl, ngl))
    
    @info ""
    @info " 3D"
    @info size(a.x)

    return a, flux
end


function myfun!(params)
    for i=1:params.nelem
        for j=1:params.ngl
            for k=1:params.ngl
                for l=1:params.ngl
                    params.a.x[i,j,k,l] = 25.0
                    params.flux.F[i,j,k,l] = 25.0
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

function fun_call(myfun!, params)
    for it=1:100
        myfun!(params)
    end 
end


function fun_call(myfun!, x, nelem, ngl)
    for it=1:100
        myfun!(x, nelem, ngl)
    end 
end
    

function main()

    nelem = 10
    ngl   =  5
    
    a, flux = allocate_operator_arrays(NSD_3D(), nelem, ngl, Float64; neqs=1)
    
    @info ""
    @info " main 3D"
    @info size(a.x)

    params = (a=a, flux=flux, nelem=nelem, ngl=ngl)

    @info " Passing single arguments"
    @time fun_call(myfun!, a.x, nelem, ngl)
    
    @info " Passing tuple"
    @time fun_call(myfun!, params)

    
end

main()
