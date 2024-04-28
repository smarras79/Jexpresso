Base.@kwdef mutable struct St_a{TFloat <: AbstractFloat}
    x::Array{TFloat} = zeros(TFloat, 0)
end


function allocate_arrays(nelem, ngl, TFloat; neqs=1)
    
    a = St_a{TFloat}(x = zeros(nelem, ngl, ngl, ngl))

    return a
end

#PASS TUPLE
function myfun!(params, lsum)

    if lsum == true
        for i=1:params.nelem
            for j=1:params.ngl
                for k=1:params.ngl
                    for l=1:params.ngl

                        params.a.x[i,j,k,l] += 25.0
                        
                    end
                end
            end
        end
    else
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
        
end

#PASS SINGLE ARGUMENTS
function myfun!(x, nelem, ngl, lsum)

    if lsum == true
        for i=1:nelem
            for j=1:ngl
                for k=1:ngl
                    for l=1:ngl
                        x[i,j,k,l] += 25.0
                    end
                end
            end
        end
    else
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
end

function fun_call(myfun!, params, lsum)
    #
    # passing TUPLE to fun_call:
    #
    for it=1:1000
        myfun!(params, lsum)
    end
    
end


function fun_call(myfun!, x, nelem, ngl, lsum)
    #
    # passing SINGLE ARGUMENTS to fun_call:
    #
    for it=1:1000
        myfun!(x, nelem, ngl, lsum)
    end
end
    

function main()

    nelem = 10
    ngl   =  5
    
    a = allocate_arrays(nelem, ngl, Float64; neqs=1)
    
    params = (a=a, nelem=nelem, ngl=ngl)
    
    @info "--------------------------------------------"
    @info " passing TUPLE vs passing SINGLE ARGUMENTS:"
    @info "--------------------------------------------"
    @info ""
    @info "--------------------------------------------"
    @info " NO SUM:"
    @info "--------------------------------------------"
    lsum = false

    @info " Passing SINGLE ARGUMENTS"
    @time fun_call(myfun!, a.x, nelem, ngl, lsum)
    @info "*********************************************"
    @info " Passing TUPLE"
    @time fun_call(myfun!, params, lsum)

    @info ""
    @info ""
    @info "--------------------------------------------"
    @info " WITH SUM:"
    @info "--------------------------------------------"
    lsum = true
    
    @info " Passing SINGLE ARGUMENTS"
    @time fun_call(myfun!, a.x, nelem, ngl, lsum)
    @info "*********************************************"
    @info " Passing TUPLE"
    @time fun_call(myfun!, params, lsum)
   
    
    
end

main()
