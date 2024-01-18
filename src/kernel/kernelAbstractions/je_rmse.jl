using Metal
using KernelAbstractions

@kernel function ka_rmse_kernel(C, A, B)
    i = @index(Global)

    a = A[i]
    b = B[i]
    KernelAbstractions.@atomic C[] += (a-b)^2
end

A = rand(Float32, 512, 512)
B = rand(Float32, 512, 512)

backend = CPU()
dA = KernelAbstractions.allocate(backend, eltype(A), size(A))
KernelAbstractions.copyto!(backend, dA, A)
dB = KernelAbstractions.allocate(backend, eltype(B), size(B))
KernelAbstractions.copyto!(backend, dB, B)
dC = KernelAbstractions.zeros(backend, eltype(A), 1)
@show typeof(dC)

k = ka_rmse_kernel(backend)
k(dC, dA, dB; ndrange=size(A))

sqrt(Array(dC)[] / length(A))
     

backend = MetalBackend()
dA = KernelAbstractions.allocate(backend, eltype(A), size(A))
KernelAbstractions.copyto!(backend, dA, A)
dB = KernelAbstractions.allocate(backend, eltype(B), size(B))
KernelAbstractions.copyto!(backend, dB, B)
dC = KernelAbstractions.zeros(backend, eltype(A), 1)
@show typeof(dC)

k = ka_rmse_kernel(backend)
k(dC, dA, dB; ndrange=size(A))

sqrt(Array(dC)[] / length(A))
