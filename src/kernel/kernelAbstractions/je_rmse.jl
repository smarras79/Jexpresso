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


A = rand(Float32, 20048, 20048)
B = rand(Float32, 20048, 20048)

backend = CPU()
dA = KernelAbstractions.allocate(backend, eltype(A), size(A))
KernelAbstractions.copyto!(backend, dA, A)
dB = KernelAbstractions.allocate(backend, eltype(B), size(B))
KernelAbstractions.copyto!(backend, dB, B)
dC = KernelAbstractions.zeros(backend, eltype(A), 1)
@show typeof(dC)

@btime ka_rmse_kernel(backend)

#k(dC, dA, dB; ndrange=size(A))
#sqrt(Array(dC)[] / length(A))
     

backend = MetalBackend()
dA = KernelAbstractions.allocate(backend, eltype(A), size(A))
KernelAbstractions.copyto!(backend, dA, A)
dB = KernelAbstractions.allocate(backend, eltype(B), size(B))
KernelAbstractions.copyto!(backend, dB, B)
dC = KernelAbstractions.zeros(backend, eltype(A), 1)

@benchmark ka_rmse_kernel(backend)

@benchmark ka_rmse_kernel_optimal1(backend)

@benchmark ka_rmse_kernel_optimal2(backend)

#k(dC, dA, dB; ndrange=size(A))
#sqrt(Array(dC)[] / length(A))
