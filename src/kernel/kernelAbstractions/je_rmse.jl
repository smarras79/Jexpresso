using Metal
using KernelAbstractions
using BenchmarkTools

@kernel function ka_rmse_kernel(C, A, B)
    i = @index(Global)

    a = A[i]
    b = B[i]
    KernelAbstractions.@atomic C[] += (a-b)^2 #atomic avoids race condition. Notice that @atomic is slow. write the kernel faster without using @atomic.
end

@kernel function ka_rmse_kernel_optimal1(C, A, B)

    #This is similar to above, but it uses `stride` to minimize @atomic operations
    
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x

    # grid-stride loop
    while i <= length(A)
        a = A[i]
        b = B[i]
        KernelAbstractions.@atomic C[] += (a-b)^2

        i += stride
    end
end

@kernel function ka_rmse_kernel_optimal2(C, A, B)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x

    # grid-stride loop
    val = zero(T)
    while i <= length(A)
        a = A[i]
        b = B[i]
        val += (a-b)^2                 # <-- changed

        i += stride
    end
    KernelAbstractions.@atomic C[] += val            # <-- changed
end

@kernel function ka_rmse_kernel_optimal3(C, A, B)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x

    # grid-stride loop to process elements on a single thread
    val = zero(T)
    while i <= length(A)
        a = A[i]
        b = B[i]
        val += (a-b)^2

        i += stride
    end

    # initialize shared memory
    thread = threadIdx().x
    threads = blockDim().x
    shared = CuDynamicSharedArray(eltype(C), (threads,))
    shared[thread] = val

    # perform a parallel reduction
    d = 1
    while d < threads
        sync_threads()
        index = 2 * d * (thread-1) + 1
        if index <= threads
            other_val = if index + d <= threads
                shared[index+d]
            else
                zero(eltype(C))
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

#=let kernel = @kernel launch=false rmse_kernel(C, A, B)
    # we need to inform the occupancy API about the shared memory usage!
    shmem(threads) = threads * sizeof(T)

    config = KernelAbstractions.launch_configuration(kernel.fun; shmem)
    threads = min(length(A), config.threads)
    blocks = cld(length(A), threads)
    blocks = min(blocks, config.blocks)

    kernel(C, A, B; threads, blocks, shmem=shmem(threads))
end
=#

A = rand(Float32, 2048, 2048)
B = rand(Float32, 2048, 2048)

backend = CPU()
dA = KernelAbstractions.allocate(backend, eltype(A), size(A))
KernelAbstractions.copyto!(backend, dA, A)
dB = KernelAbstractions.allocate(backend, eltype(B), size(B))
KernelAbstractions.copyto!(backend, dB, B)
dC = KernelAbstractions.zeros(backend, eltype(A), 1)
@show typeof(dC)

k= ka_rmse_kernel(backend)
@btime k(dC, dA, dB; ndrange=size(A))
#sqrt(Array(dC)[] / length(A))
     

backend = MetalBackend()
dA = KernelAbstractions.allocate(backend, eltype(A), size(A))
KernelAbstractions.copyto!(backend, dA, A)
dB = KernelAbstractions.allocate(backend, eltype(B), size(B))
KernelAbstractions.copyto!(backend, dB, B)
dC = KernelAbstractions.zeros(backend, eltype(A), 1)

k = ka_rmse_kernel(backend)
@btime k(dC, dA, dB; ndrange=size(A))

k = ka_rmse_kernel_optimal1(backend)
@btime k(dC, dA, dB; ndrange=size(A))

k = ka_rmse_kernel_optimal2(backend)
@btime k(dC, dA, dB; ndrange=size(A))

k = ka_rmse_kernel_optimal3(backend)
@btime k(dC, dA, dB; ndrange=size(A))
#sqrt(Array(dC)[] / length(A))
