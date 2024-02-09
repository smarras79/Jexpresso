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

@kernel function _build_rhs_gpu_2D_v0!(RHS, u, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, M, ngl)
    s = Int32(@groupsize()[1])
    #n = div(@ndrange()[1],s)#div(length(A),s)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    ip = connijk[ie,i_x,i_y]

    KernelAbstractions.@atomic RHS[ip] = Float32(0.0)
    DIM = @uniform @groupsize()[1]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM+1,DIM+1)
    #F[i_x,i_y] = Float32(1.0)*u[ip] #user_flux(u[ip])
    G = @localmem eltype(RHS) (DIM+1,DIM+1)
    #G[i_x,i_y] = Float32(1.0)*u[ip]
    #F[i_x,i_y] = flux[1]
    #G[i_x,i_y] = flux[2]
    #user_flux!(F[i_x,i_y],G[i_x,i_y],u[ip])
    #F[i_x,i_y], G[i_x,i_y] = user_flux(u[ip])
    flux = user_flux(u[ip])
    F[i_x,i_y] = flux[1]
    G[i_x,i_y] = flux[2]
    @synchronize()
    ### do numerical integration
    dFdξ = zero(Float32)
    dFdη = zero(Float32)
    dGdξ = zero(Float32)
    dGdη = zero(Float32)

    for k=1:ngl
        dFdξ += dψ[k,i_x]*F[k,i_y]
        dFdη += dψ[k,i_y]*F[i_x,k]
        dGdξ += dψ[k,i_x]*G[k,i_y]
        dGdη += dψ[k,i_y]*G[i_x,k]
    end

    dξdx_ij = dξdx[ie,i_x,i_y]
    dξdy_ij = dξdy[ie,i_x,i_y]
    dηdx_ij = dηdx[ie,i_x,i_y]
    dηdy_ij = dηdy[ie,i_x,i_y]

    dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
    dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
    KernelAbstractions.@atomic RHS[ip] -= ω[i_x]*ω[i_y]*Je[ie,i_x,i_y]*(dFdx + dGdy)/ M[ip]
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
