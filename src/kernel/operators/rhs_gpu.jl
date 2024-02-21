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

@kernel function _build_rhs_gpu_2D_v0!(RHS, u, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, M, ngl, neq)
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
    G = @localmem eltype(RHS) (DIM+1,DIM+1)
    S = @localmem eltype(RHS) (DIM+1,DIM+1)

    flux = user_flux(@view(u[ip,1:neq]))
    
    #source = user_source(u[ip,:],x[ip],y[ip])
    @synchronize()
    ### do numerical integration
    for ieq =1:neq
        F[i_x,i_y] = flux[ieq]
        G[i_x,i_y] = flux[neq+ieq]
        #S[i_x,i_y] = source[ieq]
        @synchronize()
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
    KernelAbstractions.@atomic RHS[ip,ieq] -= ω[i_x]*ω[i_y]*Je[ie,i_x,i_y]*((dFdx + dGdy))/M[ip]#- S[i_x,i_y])/ M[ip]
    end
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
