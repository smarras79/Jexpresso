@kernel function _build_rhs_gpu_v0!(RHS, u, connijk, dψ, ω, M, ngl)
    s = Int32(@groupsize()[1])
    #n = div(@ndrange()[1],s)#div(length(A),s)
    ig = @index(Group, Linear)
    il = @index(Local, Linear)
    ip = connijk[ig,il,1]
    
    KernelAbstractions.@atomic RHS[ip] = Float32(0.0)
    DIM = @uniform @groupsize()[1]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM+1,1)
    F[il] = Float32(1.0)*u[ip] #user_flux(u[ip])
    @synchronize()
    ### do numerical integration
    dFdxi = zero(Float32)
    for k=1:ngl
        dFdxi += dψ[k,il]*F[k]
    end
    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
    KernelAbstractions.@atomic RHS[ip] -= ω[il]*dFdxi/ M[ip]
end

#=backend = MetalBackend()
ngl = 5
ne = 100
ω = KernelAbstractions.ones(backend, Float32, ngl)
dψ = KernelAbstractions.ones(backend, Float32, ngl,ngl)
connijk = KernelAbstractions.ones(backend, Int32, ne,ngl,1)
M = KernelAbstractions.ones(backend, Float32, ne*ngl)
u = KernelAbstractions.ones(backend, Float32, ne*ngl)
RHS = KernelAbstractions.ones(backend, Float32, ne*ngl)
k = _build_rhs_gpu_v0!(backend,(ngl))
for i=1:ne*ngl
    connijk[i] = ne*ngl+1 -i
end



k(RHS, u, connijk , dψ, ω, M, ngl; ndrange = ne*ngl,workgroupsize = ngl)=#
