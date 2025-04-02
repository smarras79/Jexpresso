@kernel function ka_rmse_kernel(C, A, B)
    i = @index(Global)

    a = A[i]
    b = B[i]
    KernelAbstractions.@atomic C[] += (a-b)^2 #atomic avoids race condition. Notice that @atomic is slow. write the kernel faster without using @atomic.
end

@kernel function ka_rmse_kernel_optimal1(C, A, B)

    #This is similar to above, but it uses `stride` to minimize @atomic operations   
    s = Int32(@groupsize()[1])
    n = div(@ndrange()[1]*@ndrange()[2],s)
    ig = @index(Group, Linear)
    il = @index(Local, Linear)
    i=(ig-1)*s + il
    stride = s*n
    while i <= length(A)
        a = A[i]
        b = B[i]
        KernelAbstractions.@atomic C[] += (a-b)^2

        i += stride
    end
    
end

@kernel function ka_rmse_kernel_optimal2(C, A, B)
    #i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    #stride = gridDim().x * blockDim().x

    s = Int32(@groupsize()[1])
    n = div(@ndrange()[1]*@ndrange()[2],s)
    ig = @index(Group, Linear)
    il = @index(Local, Linear)
    i=(ig-1)*s + il
    stride = s*n

    # grid-stride loop
    val = zero(Float32)
    while i <= length(A)
        a = A[i]
        b = B[i]
        val += (a-b)^2                 # <-- changed

        i += stride
    end
    KernelAbstractions.@atomic C[] += val            # <-- changed
end

@kernel function ka_rmse_kernel_optimal3(C, A, B)
    s = Int32(@groupsize()[1])
    n = div(@ndrange()[1],s)#div(length(A),s)
    ig = @index(Group, Linear)
    il = @index(Local, Linear)
    i=(ig-1)*s + il
    stride = s*n
    
    # grid-stride loop to process elements on a single thread
    val = zero(Float32)
    while i <= length(A)
        a = A[i]
        b = B[i]
        val += (a-b)^2

        i += stride
    end

    # initialize shared memory
    thread = il
    threads = s
    DIM = @uniform @groupsize()[1]
    shared = @localmem eltype(C) (DIM+1,1) 
     
    shared[thread] = val
    #@info shared[thread], val
    # perform a parallel reduction
    d = Int32(1)
    #@info shared
    while d < threads
        
        @synchronize()
        index = 2 * d * (thread-1) + 1
        if index <= threads
            other_val = zero(eltype(C))
            if (index + d <= threads)
                other_val = shared[index+d]
            end
            shared[index] = shared[index] + other_val
        end
        d *= 2
    end

    # load the final value on the first thread
    if thread == 1
        val = shared[thread]
        KernelAbstractions.@atomic C[] += val
    end
end


@kernel function prototype_rhs!(rhs,omega,dpsi,x,M,u,ngl)
    s = Int32(@groupsize()[1])
    n = div(@ndrange()[1],s)#div(length(A),s)
    ig = @index(Group, Linear)
    il = @index(Local, Linear)
    i=(ig-1)*s + il
    DIM = @uniform @groupsize()[1]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(rhs) (DIM+1,1)
    F[il] = user_flux(u[i])
    @synchronize()
    ### do numerical integration
    dFdxi = zero(Float32)
    for k=1:ngl
       dFdxi += dpsi[k,il]*F[k]
    end
    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step 
    KernelAbstractions.@atomic rhs[i] -= omega[il]*dFdxi/ M[i]
end

function user_flux(u)
  return 2.0*u
end


N=512
A = rand(Float32, N, N)
B = rand(Float32, N, N)

backend = CPU()
dA = KernelAbstractions.allocate(backend, eltype(A), size(A))
KernelAbstractions.copyto!(backend, dA, A)
dB = KernelAbstractions.allocate(backend, eltype(B), size(B))
KernelAbstractions.copyto!(backend, dB, B)
dC = KernelAbstractions.zeros(backend, eltype(A), 1)
@show typeof(dC)

@info "CPU"
k= ka_rmse_kernel(backend)
@time k(dC, dA, dB; ndrange=size(A))
#sqrt(Array(dC)[] / length(A))
KernelAbstractions.synchronize(backend)     
@info dC
#dC .= 0.0
@info "GPU 1"

backend = MetalBackend()
#backend = CPU()
dA = KernelAbstractions.allocate(backend, eltype(A), size(A))
KernelAbstractions.copyto!(backend, dA, A)
dB = KernelAbstractions.allocate(backend, eltype(B), size(B))
KernelAbstractions.copyto!(backend, dB, B)
dC = KernelAbstractions.zeros(backend, eltype(A), 1)

k = ka_rmse_kernel(backend)
@time k(dC, dA, dB; ndrange=size(A),workgroupsize = 512)
KernelAbstractions.synchronize(backend)

@info dC
dC .= Float32(0.0)

@info "GPU 2"
k = ka_rmse_kernel_optimal1(backend)
@time k(dC, dA, dB; ndrange=size(A),workgroupsize = 512)
KernelAbstractions.synchronize(backend)
@info dC
dC .= Float32(0.0)

@info "GPU 3"
k = ka_rmse_kernel_optimal2(backend)
@time k(dC, dA, dB; ndrange=size(A),workgroupsize = 512)
KernelAbstractions.synchronize(backend)
@info dC
dC .= Float32(0.0)


@info "GPU 4"
k = ka_rmse_kernel_optimal3(backend,(256))
@time k(dC, dA, dB; ndrange=256,workgroupsize = 256)
KernelAbstractions.synchronize(backend)
@info dC

ngl = 5
ne = 100
omega = KernelAbstractions.ones(backend, eltype(A), ngl)
dpsi = KernelAbstractions.ones(backend, eltype(A), ngl,ngl)
x = KernelAbstractions.ones(backend, eltype(A), ne,ngl)
M = KernelAbstractions.ones(backend, eltype(A), ne*ngl)
u = KernelAbstractions.ones(backend, eltype(A), ne*ngl)
rhs = KernelAbstractions.ones(backend, eltype(A), ne*ngl)
k = prototype_rhs!(backend,(ngl))
@time k(rhs,omega,dpsi,x,M,u,ngl; ndrange = ne*ngl,workgroupsize = ngl)
