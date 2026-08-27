using Quadmath

#---------------------------------------------------------------------------------------
# Dummy Laguerre arguments.
#
# These used to be `zeros(TFloat,1,1,1)` written directly in the keyword-argument
# default of every `filter!` method.  Two problems with that:
#
#   1. the default expression is evaluated on EVERY call that does not pass the
#      Laguerre arrays (i.e. every call of every non-Laguerre run), so each RHS
#      evaluation heap-allocated two throw-away arrays that are never read;
#   2. `TFloat` is a *non-const* global (see src/Jexpresso.jl and the
#      `global TFloat = Float32` in src/io/mod_inputs.jl), so `zeros(TFloat,...)`
#      is type-unstable: the keyword arguments come out as `Any`, which forces
#      boxing and a dynamic dispatch on entry to `filter!`.
#
# Hoisting them to `const` globals removes both the allocation and the
# instability.  They are zero-sized because they are never indexed: they only
# exist so the Laguerre branch type-checks.
#---------------------------------------------------------------------------------------
const FILTER_NO_CONN_LAG_2D = zeros(Int64,   0, 0, 0)
const FILTER_NO_JE_LAG_2D   = zeros(Float64, 0, 0, 0)
const FILTER_NO_CONN_LAG_3D = zeros(Int64,   0, 0, 0, 0)
const FILTER_NO_JE_LAG_3D   = zeros(Float64, 0, 0, 0, 0)


#---------------------------------------------------------------------------------------
# NOTE on the element-local scratch array `b` (params.b):
#
# `b` is dimensioned (nelem, ngl^d, neqs) — for a 3D LES mesh that is tens to
# hundreds of MB.  It used to be
#
#     fill!(b, 0)                       # in reset_filters!, every RHS evaluation
#     b[e,i,j,k,m] += fqf[...]*w*Je     # here
#     DSS_rhs!(B, b, ...)               # immediately reads it back and scatters into B
#
# i.e. the array was memset, written once, and read once — and `e`, the *first*
# (fastest varying) index, is the *outermost* loop, so every scalar store lands
# a full `nelem*sizeof(T)` stride away from the previous one.  On a 7680-element
# mesh that is a fresh cache line (and often a fresh page) per Float64.
#
# `b` is only genuinely needed when `ladapt == true`, because `DSS_nc_gather_rhs!`
# consumes it.  For every other run the direct-stiffness sum can simply be done
# in place inside the element loop, exactly the way `filter_gpu_2d!`/`filter_gpu_3d!`
# already do it on the GPU.  That removes the memset, the array traffic and the
# separate `DSS_rhs!` sweep.
#
# The scatter order is unchanged (elements are still visited in increasing `e`
# for each (point, equation) pair), so the assembled result is bit-for-bit
# identical to the previous implementation.
#
# When `b` *is* needed, it is now written with `=` instead of `+=`; every entry
# is assigned exactly once per call, so `fill!(b, 0)` was redundant there too
# and has been dropped from `reset_filters!`.
#---------------------------------------------------------------------------------------

function filter!(u, params, t, uaux, connijk, Je, SD::NSD_2D, ::TOTAL;
                 connijk_lag = FILTER_NO_CONN_LAG_2D, Je_lag = FILTER_NO_JE_LAG_2D, ladapt = false)

    # NOTE the Int64(...) conversions: St_mesh declares npoin/nelem/ngl as
    # ::Union{TInt, Missing}, so reading them straight into a loop bound leaves
    # every trip count union-typed and forces a union-split on each loop entry.
    # The old code re-read `params.mesh.ngl` inside the innermost tensor
    # contractions, paying for that on every single iteration. Hoisting to
    # concrete Int64 locals is the same typed-barrier trick _sphere_filter_kernel!
    # already uses.
    npoin = Int64(params.mesh.npoin)
    nelem = Int64(params.mesh.nelem)
    ngl   = Int64(params.mesh.ngl)
    neqs  = Int64(params.neqs)
    ω     = params.ω
    qe    = params.qp.qe

    q_t   = params.q_t
    q_ti  = params.q_ti
    fx    = params.fx
    fy_t  = params.fy_t
    fqf   = params.fqf
    b     = params.b
    B     = params.B
    AD    = params.AD

    TF    = eltype(B)

    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)

    ## Subtract reference state from ALL prognostic variables (ρ, ρu, ρv, ρE
    ## for TOTAL) before filtering, so the modal filter acts only on the
    ## perturbation from the IC. Otherwise the filter slowly smooths the
    ## reference state itself and leaks ringing near sharp ICs (e.g. KH shear
    ## layer in ρE). This now matches the behaviour of filter_gpu_2d!.
    @views uaux[:,1:neqs] .-= qe[:,1:neqs]

    ## Loop through the elements
    @inbounds for e = 1:nelem

        for j = 1:ngl, i = 1:ngl
            ip = connijk[e,i,j]
            for m = 1:neqs
                q_t[m,i,j] = uaux[ip,m]
            end
        end

        ### Construct local derivatives for prognostic variables
        for m = 1:neqs

            ## KSI Derivative
            for j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for k = 1:ngl
                    acc += fx[i,k] * q_t[m,k,j]
                end
                q_ti[i,j] = acc
            end

            ## ETA Derivative
            for j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for k = 1:ngl
                    acc += q_ti[i,k] * fy_t[k,j]
                end
                fqf[m,i,j] = acc
            end
        end

        ## Do Numerical Integration (+ DSS fused in, unless AMR needs `b`)
        if ladapt
            for j = 1:ngl, i = 1:ngl
                for m = 1:neqs
                    # NOTE: original left-to-right association kept on purpose,
                    # see the 3D comment -- hoisting ω[i]*ω[j]*Je changes the
                    # result at the ulp level.
                    b[e,i,j,m] = fqf[m,i,j] * ω[i]*ω[j]*Je[e,i,j]
                end
            end
        else
            for j = 1:ngl, i = 1:ngl
                ip = connijk[e,i,j]
                for m = 1:neqs
                    B[ip,m] += fqf[m,i,j] * ω[i]*ω[j]*Je[e,i,j]
                end
            end
        end
    end

    if ladapt == true
        DSS_nc_gather_rhs!(B, SD, params.QT, b, connijk, params.mesh.poin_in_edge,
                           params.mesh.non_conforming_facets, params.mesh.cip, params.mesh.pip, params.mesh.lfid, params.mesh.half1, params.mesh.half2,
                           params.mesh.non_conforming_facets_parents_ghost, params.mesh.cip_pg, params.mesh.lfid_pg, params.mesh.half1_pg, params.mesh.half2_pg,
                           params.q_el, params.q_el_pro, params.L_1, params.L_2, params.q_ghost_p,
                           params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPc_list_pg,
                           params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost, params.mesh.pgip_owner, params.mesh.pgip_local,
                           ngl-1, neqs, params.interp)
        DSS_rhs!(B, b, connijk, nelem, ngl, neqs, SD, AD)
    end

    DSS_global_RHS!(B, params.g_dss_cache, neqs)
    for ieq = 1:neqs
        divide_by_mass_matrix!(@view(B[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
    end

    @views uaux[:,1:neqs] .= B[:,1:neqs]
    @views uaux[:,1:neqs] .+= qe[:,1:neqs]

    uaux2u!(u, @view(uaux[:,:]), neqs, npoin)
end


function filter!(u, params, t, uaux, connijk, Je, SD::NSD_2D, ::PERT;
                 connijk_lag = FILTER_NO_CONN_LAG_2D, Je_lag = FILTER_NO_JE_LAG_2D, ladapt = false)

    # NOTE the Int64(...) conversions: St_mesh declares npoin/nelem/ngl as
    # ::Union{TInt, Missing}, so reading them straight into a loop bound leaves
    # every trip count union-typed and forces a union-split on each loop entry.
    # The old code re-read `params.mesh.ngl` inside the innermost tensor
    # contractions, paying for that on every single iteration. Hoisting to
    # concrete Int64 locals is the same typed-barrier trick _sphere_filter_kernel!
    # already uses.
    npoin = Int64(params.mesh.npoin)
    nelem = Int64(params.mesh.nelem)
    ngl   = Int64(params.mesh.ngl)
    neqs  = Int64(params.neqs)
    ω     = params.ω

    q_t   = params.q_t
    q_ti  = params.q_ti
    fx    = params.fx
    fy_t  = params.fy_t
    fqf   = params.fqf
    b     = params.b
    B     = params.B
    AD    = params.AD

    TF    = eltype(B)

    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)

    ## Loop through the elements
    @inbounds for e = 1:nelem

        for j = 1:ngl, i = 1:ngl
            ip = connijk[e,i,j]
            for m = 1:neqs
                q_t[m,i,j] = uaux[ip,m]
            end
        end

        ### Construct local derivatives for prognostic variables
        for m = 1:neqs

            ## KSI Derivative
            for j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for k = 1:ngl
                    acc += fx[i,k] * q_t[m,k,j]
                end
                q_ti[i,j] = acc
            end

            ## ETA Derivative
            for j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for k = 1:ngl
                    acc += q_ti[i,k] * fy_t[k,j]
                end
                fqf[m,i,j] = acc
            end
        end

        ## Do Numerical Integration (+ DSS fused in, unless AMR needs `b`)
        if ladapt
            for j = 1:ngl, i = 1:ngl
                for m = 1:neqs
                    # NOTE: original left-to-right association kept on purpose,
                    # see the 3D comment -- hoisting ω[i]*ω[j]*Je changes the
                    # result at the ulp level.
                    b[e,i,j,m] = fqf[m,i,j] * ω[i]*ω[j]*Je[e,i,j]
                end
            end
        else
            for j = 1:ngl, i = 1:ngl
                ip = connijk[e,i,j]
                for m = 1:neqs
                    B[ip,m] += fqf[m,i,j] * ω[i]*ω[j]*Je[e,i,j]
                end
            end
        end
    end

    if ladapt
        DSS_rhs!(B, b, connijk, nelem, ngl, neqs, SD, AD)
    end

    if (params.laguerre)

        ngr        = Int64(params.mesh.ngr)
        nelem_semi = Int64(params.mesh.nelem_semi_inf)
        ω_lag      = params.ω_lag
        q_t_lag    = params.q_t_lag
        q_ti_lag   = params.q_ti_lag
        fqf_lag    = params.fqf_lag
        fy_t_lag   = params.fy_t_lag
        B_lag      = params.B_lag

        @inbounds for e = 1:nelem_semi

            for j = 1:ngr, i = 1:ngl
                ip = connijk_lag[e,i,j]
                for m = 1:neqs
                    q_t_lag[m,i,j] = uaux[ip,m]
                end
            end

            ### Construct local derivatives for prognostic variables
            for m = 1:neqs

                for j = 1:ngr, i = 1:ngl
                    acc = zero(TF)
                    for k = 1:ngl
                        acc += fx[i,k] * q_t_lag[m,k,j]
                    end
                    q_ti_lag[i,j] = acc
                end

                ## ETA Derivative
                for j = 1:ngr, i = 1:ngl
                    acc = zero(TF)
                    for k = 1:ngr
                        acc += q_ti_lag[i,k] * fy_t_lag[k,j]
                    end
                    fqf_lag[m,i,j] = acc
                end
            end

            ## Numerical integration + DSS of the Laguerre block, fused.
            ## (2D Laguerre + AMR is not supported, see the @mystop below, so
            ##  `b_lag` is never needed and is not written at all.)
            for j = 1:ngr, i = 1:ngl
                ip = connijk_lag[e,i,j]
                for m = 1:neqs
                    B_lag[ip,m] += fqf_lag[m,i,j] * ω[i]*ω_lag[j]*Je_lag[e,i,j]
                end
            end
        end

        if ladapt == true
            @mystop("not yet implemented for 2D Laguerre with AMR!")
        end

        B .+= B_lag
    end

    DSS_global_RHS!(B, params.g_dss_cache, neqs)

    for ieq = 1:neqs
        divide_by_mass_matrix!(@view(B[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
    end
    @views uaux[:,1:neqs] .= B[:,1:neqs]

    uaux2u!(u, @view(uaux[:,:]), neqs, npoin)
end


function filter!(u, params, t, uaux, connijk, Je, SD::NSD_3D, ::PERT;
                 connijk_lag = FILTER_NO_CONN_LAG_3D, Je_lag = FILTER_NO_JE_LAG_3D, ladapt = false)

    # NOTE the Int64(...) conversions: St_mesh declares npoin/nelem/ngl as
    # ::Union{TInt, Missing}, so reading them straight into a loop bound leaves
    # every trip count union-typed and forces a union-split on each loop entry.
    # The old code re-read `params.mesh.ngl` inside the innermost tensor
    # contractions, paying for that on every single iteration. Hoisting to
    # concrete Int64 locals is the same typed-barrier trick _sphere_filter_kernel!
    # already uses.
    npoin = Int64(params.mesh.npoin)
    nelem = Int64(params.mesh.nelem)
    ngl   = Int64(params.mesh.ngl)
    neqs  = Int64(params.neqs)
    ω     = params.ω

    q_t   = params.q_t
    q_ti  = params.q_ti
    fx    = params.fx
    q_tij = params.q_tij
    fy_t  = params.fy_t
    fqf   = params.fqf
    fz_t  = params.fz_t
    b     = params.b
    B     = params.B

    QT    = params.QT
    AD    = params.AD

    TF    = eltype(B)

    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)

    ## Loop through the elements
    @inbounds for e = 1:nelem

        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = connijk[e,i,j,k]
            for m = 1:neqs
                q_t[m,i,j,k] = uaux[ip,m]
            end
        end

        ### Construct local derivatives for prognostic variables
        for m = 1:neqs

            ## KSI derivative
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for l = 1:ngl
                    acc += fx[i,l] * q_t[m,l,j,k]
                end
                q_ti[i,j,k] = acc
            end

            ## ETA derivative
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for l = 1:ngl
                    acc += q_ti[i,l,k] * fy_t[l,j]
                end
                q_tij[i,j,k] = acc
            end

            ## ZETA derivative
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for l = 1:ngl
                    acc += q_tij[i,j,l] * fz_t[l,k]
                end
                fqf[m,i,j,k] = acc
            end
        end

        ## Do Numerical Integration (+ DSS fused in, unless AMR needs `b`)
        if ladapt
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                for m = 1:neqs
                    # NOTE: the product is written out in the original
                    # left-to-right order (fqf first) on purpose -- hoisting
                    # ω[i]*ω[j]*ω[k]*Je out of the m loop reassociates the
                    # multiplications and perturbs the result at the ulp level.
                    b[e,i,j,k,m] = fqf[m,i,j,k] * ω[i]*ω[j]*ω[k]*Je[e,i,j,k]
                end
            end
        else
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip = connijk[e,i,j,k]
                for m = 1:neqs
                    B[ip,m] += fqf[m,i,j,k] * ω[i]*ω[j]*ω[k]*Je[e,i,j,k]
                end
            end
        end
    end

    if ladapt == true
        DSS_nc_gather_rhs!(B, SD, QT, b,
                           params.mesh.non_conforming_facets,
                           params.mesh.non_conforming_facets_parents_ghost, params.cache_ghost_p,
                           params.q_el, params.q_el_pro, params.q_ghost_p,
                           params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPc_list_pg,
                           params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost,
                           params.mesh.pgip_local,
                           ngl-1, neqs, params.interp)
        DSS_rhs!(B, b, connijk, nelem, ngl, neqs, SD, AD)
    end

    DSS_global_RHS!(B, params.g_dss_cache, neqs)
    for ieq = 1:neqs
        divide_by_mass_matrix!(@view(B[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
        if ladapt == true
            DSS_nc_scatter_rhs!(@view(B[:,ieq]), SD, QT,
                                params.mesh.non_conforming_facets,
                                params.mesh.non_conforming_facets_children_ghost, params.cache_ghost_c,
                                params.q_el, params.q_el_pro, params.q_ghost_c,
                                params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPp_list_cg,
                                params.mesh.gip2ip, params.mesh.cgip_local,
                                ngl-1, params.interp)
        end
    end

    @views uaux[:,1:neqs] .= B[:,1:neqs]

    uaux2u!(u, @view(uaux[:,:]), neqs, npoin)
end


function filter!(u, params, t, uaux, connijk, Je, SD::NSD_3D, ::TOTAL;
                 connijk_lag = FILTER_NO_CONN_LAG_3D, Je_lag = FILTER_NO_JE_LAG_3D, ladapt = false)

    # NOTE the Int64(...) conversions: St_mesh declares npoin/nelem/ngl as
    # ::Union{TInt, Missing}, so reading them straight into a loop bound leaves
    # every trip count union-typed and forces a union-split on each loop entry.
    # The old code re-read `params.mesh.ngl` inside the innermost tensor
    # contractions, paying for that on every single iteration. Hoisting to
    # concrete Int64 locals is the same typed-barrier trick _sphere_filter_kernel!
    # already uses.
    npoin = Int64(params.mesh.npoin)
    nelem = Int64(params.mesh.nelem)
    ngl   = Int64(params.mesh.ngl)
    neqs  = Int64(params.neqs)
    ω     = params.ω

    qe    = params.qp.qe

    q_t   = params.q_t
    q_ti  = params.q_ti
    fx    = params.fx
    q_tij = params.q_tij
    fy_t  = params.fy_t
    fqf   = params.fqf
    fz_t  = params.fz_t
    b     = params.b
    B     = params.B

    QT    = params.QT
    AD    = params.AD

    TF    = eltype(B)

    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)

    ## Subtract reference state from ALL prognostic variables (ρ, ρu, ρv, ρw,
    ## ρE for TOTAL) before filtering — see 2D ::TOTAL comment for rationale.
    @views uaux[:,1:neqs] .-= qe[:,1:neqs]

    ## Loop through the elements
    @inbounds for e = 1:nelem

        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = connijk[e,i,j,k]
            for m = 1:neqs
                q_t[m,i,j,k] = uaux[ip,m]
            end
        end

        ### Construct local derivatives for prognostic variables
        for m = 1:neqs

            ## KSI derivative
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for l = 1:ngl
                    acc += fx[i,l] * q_t[m,l,j,k]
                end
                q_ti[i,j,k] = acc
            end

            ## ETA derivative
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for l = 1:ngl
                    acc += q_ti[i,l,k] * fy_t[l,j]
                end
                q_tij[i,j,k] = acc
            end

            ## ZETA derivative
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                acc = zero(TF)
                for l = 1:ngl
                    acc += q_tij[i,j,l] * fz_t[l,k]
                end
                fqf[m,i,j,k] = acc
            end
        end

        ## Do Numerical Integration (+ DSS fused in, unless AMR needs `b`)
        if ladapt
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                for m = 1:neqs
                    # NOTE: the product is written out in the original
                    # left-to-right order (fqf first) on purpose -- hoisting
                    # ω[i]*ω[j]*ω[k]*Je out of the m loop reassociates the
                    # multiplications and perturbs the result at the ulp level.
                    b[e,i,j,k,m] = fqf[m,i,j,k] * ω[i]*ω[j]*ω[k]*Je[e,i,j,k]
                end
            end
        else
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip = connijk[e,i,j,k]
                for m = 1:neqs
                    B[ip,m] += fqf[m,i,j,k] * ω[i]*ω[j]*ω[k]*Je[e,i,j,k]
                end
            end
        end
    end

    if ladapt == true
        DSS_nc_gather_rhs!(B, SD, QT, b,
                           params.mesh.non_conforming_facets,
                           params.mesh.non_conforming_facets_parents_ghost, params.cache_ghost_p,
                           params.q_el, params.q_el_pro, params.q_ghost_p,
                           params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPc_list_pg,
                           params.mesh.ip2gip, params.mesh.gip2ip, params.mesh.pgip_ghost,
                           params.mesh.pgip_local,
                           ngl-1, neqs, params.interp)
        DSS_rhs!(B, b, connijk, nelem, ngl, neqs, SD, AD)
    end

    DSS_global_RHS!(B, params.g_dss_cache, neqs)
    for ieq = 1:neqs
        divide_by_mass_matrix!(@view(B[:,ieq]), params.vaux, params.Minv, neqs, npoin, AD)
        if ladapt == true
            DSS_nc_scatter_rhs!(@view(B[:,ieq]), SD, QT,
                                params.mesh.non_conforming_facets,
                                params.mesh.non_conforming_facets_children_ghost, params.cache_ghost_c,
                                params.q_el, params.q_el_pro, params.q_ghost_c,
                                params.mesh.IPc_list, params.mesh.IPp_list, params.mesh.IPp_list_cg,
                                params.mesh.gip2ip, params.mesh.cgip_local,
                                ngl-1, params.interp)
        end
    end

    @views uaux[:,1:neqs] .= B[:,1:neqs]
    @views uaux[:,1:neqs] .+= qe[:,1:neqs]

    uaux2u!(u, @view(uaux[:,:]), neqs, npoin)
end




@kernel function filter_gpu_2d!(u, qe, B, fx, fy_t, Je, ω_x, ω_y, connijk, Minv, n_x, n_y, neqs,lpert)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i = il[1]
    @inbounds j = il[2]
    @inbounds ip = connijk[ie,i,j]

    #define local arrays for element based filtering
    DIM1  = @uniform @groupsize()[1]
    DIM2  = @uniform @groupsize()[2]
    q_t  = @localmem eltype(u) (DIM1+1,DIM2+1)
    q_ti = @localmem eltype(u) (DIM1+1,DIM2+1) 
    fqf  = @localmem eltype(u) (DIM1+1,DIM2+1)

    for m=1:neqs
        if (lpert)
            @inbounds q_t[i,j] = u[ip,m]
        else
            @inbounds q_t[i,j] = u[ip,m] - qe[ip,m]
        end
        @synchronize()
        q_ti[i,j] = zero(eltype(u)) 
        for k=1:n_x
            @inbounds q_ti[i,j] += fx[i,k] * q_t[k,j]
        end
        
        @synchronize()
        fqf[i,j] = zero(eltype(u))
        for k=1:n_y
            @inbounds fqf[i,j] += q_ti[i,k] * fy_t[k,j]
        end

        @inbounds KernelAbstractions.@atomic B[ip,m] += ω_x[i]*ω_y[j]*Je[ie,i,j]*fqf[i,j] * Minv[ip]
    end
end

@kernel function filter_gpu_3d!(u, qe, B, fx, fy_t, fz_t, Je, ω_x, ω_y, ω_z, connijk, Minv, n_x, n_y, n_z, neqs, lpert)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i = il[1]
    @inbounds j = il[2]
    @inbounds k = il[3]
    @inbounds ip = connijk[ie,i,j,k]

    #define local arrays for element based filtering
    DIM   = @uniform @groupsize()[1]
    q_t   = @localmem eltype(u) (DIM+1,DIM+1,DIM+1)
    q_ti  = @localmem eltype(u) (DIM+1,DIM+1,DIM+1)
    q_tij = @localmem eltype(u) (DIM+1,DIM+1,DIM+1)
    fqf   = @localmem eltype(u) (DIM+1,DIM+1,DIM+1)

    for m=1:neqs
        if (lpert)
            @inbounds q_t[i,j,k] = u[ip,m]
        else
            @inbounds q_t[i,j,k] = u[ip,m] - qe[ip,m]
        end
        @synchronize()
        q_ti[i,j,k] = zero(eltype(u))
        for l=1:n_x
            @inbounds q_ti[i,j,k] += fx[i,l] * q_t[l,j,k]
        end

        @synchronize()
        q_tij[i,j,k] = zero(eltype(u))
        for l=1:n_y
            @inbounds q_tij[i,j,k] += q_ti[i,l,k] * fy_t[l,j]
        end

        @synchronize()
        fqf[i,j,k] = zero(eltype(u))
        for l=1:n_z
            # BUGFIX: this contracted q_ti (the ksi-filtered field) instead of
            # q_tij (ksi- AND eta-filtered), so the GPU 3D filter applied the
            # eta filter and then threw it away. The CPU path has always used
            # q_tij here; this makes filter_gpu_3d! agree with filter!.
            @inbounds fqf[i,j,k] += q_tij[i,j,l] * fz_t[l,k]
        end


        @inbounds KernelAbstractions.@atomic B[ip,m] += ω_x[i]*ω_y[j]*ω_z[k]*Je[ie,i,j,k]*fqf[i,j,k]*Minv[ip]
    end 
end


function init_filter(nop,xgl,mu_x,mesh,inputs, rank)
    
    if rank == 0
        println(" # Legendre filter")
    end
    
    f = zeros(TFloat,nop+1,nop+1)
    weight = ones(Float64,nop+1)
    exp_alpha = 36
    exp_order = 64
    quad_alpha = 1.0
    quad_order = (nop+1)/3
    erf_alpha = 0.0
    erf_order = 12
    Legendre = St_Legendre{Float128}(0.0,0.0,0.0,0.0)  
    leg = zeros(Float128,nop+1,nop+1)
    ## Legendre Polynomial matrix
    if (nop+1 == mesh.ngl)
        
        for i = 1:nop+1
            ξ = xgl[i]
            for j = 1:nop+1
                jj = j - 1
                LegendreAndDerivativeAndQ!(Legendre, jj, ξ)      
                leg[i,j] = Legendre.legendre
            end
        end
        
        ### Heirarchical Modal Legendre Basis
        leg2 = zeros(Float128,nop+1,nop+1)
        leg2 .= leg
        for i=1:nop+1
            ξ = xgl[i]
            leg2[i,1] = 0.5*(1 - ξ)
            if (nop +1 > 1)
                leg2[i,2] = 0.5*(1 + ξ)
                for j=3:nop+1
                    leg2[i,j] = leg[i,j] - leg[i,j-2]
                end
            end
        end
    elseif (nop+1 == mesh.ngr)
        
        Laguerre = St_Laguerre(Polynomial(Float128(2.0)),Polynomial(Float128(2.0)),Polynomial(Float128(2.0)),Polynomial(Float128(2.0)))
        for i=1:nop+1
            ξ = xgl[i]
            for j=1:nop+1
                jj = j-1
                ScaledLaguerreAndDerivative!(jj,Laguerre,1.0,inputs[:backend])
                leg[i,j] = Laguerre.Laguerre(ξ)
            end
        end
        
        ### Scaled Laguerre Basis
        leg2 = zeros(Float128,nop+1,nop+1)
        leg2 .= leg
        for i=1:nop+1
            ξ = xgl[i]
            leg2[i,1] = exp(-ξ/2)
            if (nop +1 > 1)
                for j=2:nop+1
                    leg2[i,j] = exp(-ξ/2)*(leg[i,j] - leg[i,j-1])
                end
            end
        end
        
    end
    #### Compute Inverse Matrix
    leg_inv = zeros(Float128,nop+1,nop+1)
    leg_inv .= leg2
    ierr = 0
    gaujordf!(leg_inv,nop+1,ierr)
    if (ierr != 0)
        println(" # Error in GAUJORDF in FILTER INIT")
        @info "ierr", ierr
        exit
    end
    
    ## Compute Boyd-Vandeven (ERF-LOG) Transfer function
    filter_type = inputs[:filter_type]
    if (filter_type == "erf")   
        if rank == 0
            println(" # erf filtering on")
        end
        
        for k=1:nop+1
            # Boyd filter
            weight[k] = vandeven_modal(k,nop+1,erf_order)
        end
    elseif (filter_type == "quad")
        if rank == 0
            println(" # quadratic filtering on")
        end
        mode_filter = floor(quad_order)   
        k0 = Int64(nop+1 - mode_filter)
        xmode2 = mode_filter*mode_filter
        weight .= 1
        for k=k0+1:nop+1
            amp = quad_alpha*(k-k0)*(k-k0)/(xmode2)
            weight[k] = 1.0 - amp
        end
    elseif (filter_type == "exp")
        if rank == 0
            println(" # exponential filtering on")
        end
        for k=1:nop+1
            weight[k] = exp(-exp_alpha*(Float64(k-1)/nop)^exp_order)
        end
    end

    ### Use Laguerre Weights perhaps???
    #=if (nop +1 == mesh.ngr)
    ξω2 = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta])
    ξ2,ω2 = ξω2.ξ, ξω2.ω 
    for k=1:nop+1
    weight[k] = ω2[k]
    end
    end=# ####This doesn't do a good job 

    ## Construct 1D Filter matrix
    for i=1:nop+1
        for j=1:nop+1
            sum = 0
            for k=1:nop+1
                sum = sum + leg2[i,k] * weight[k] * leg_inv[k,j]
            end
            f[i,j] = mu_x * sum
        end
        f[i,i] = f[i,i] + (1.0 - mu_x)
    end
    return f
end


function gaujordf!(a,n,ierr)
    
    ## Initialize
    ierr = 0
    eps = 1.0e-9
    ipiv = zeros(Int64,n)
    indr = zeros(Int64,n)
    indc = zeros(Int64,n)  

    for i=1:n
        big = 0

        ## Pivot Search
        irow = -1
        icol = -1
        
        for j=1:n
            if (ipiv[j] != 1)
                for k = 1:n
                    if (ipiv[k] == 0)
                        if (abs(a[j,k]) >= big)
                            big = abs(a[j,k])
                            irow = j
                            icol = k
                        end
                    elseif (ipiv[k] > 1)
                        ierr = -ipiv
                        return nothing
                    end
                end
            end
        end
        
        ipiv[icol] = ipiv[icol] + 1

        ## Swap rows
        
        if (irow != icol)
            for l=1:n
                dum = a[irow,l]
                a[irow,l] = a[icol,l]
                a[icol,l] = dum
            end
        end
        indr[i] = irow
        indc[i] = icol
        if (abs(a[icol,icol]) < eps)
            @info "small Gauss Jordan Pivot:", icol, a[icol,icol]
            ierr = icol
            return nothing
        end
        piv = 1.0/a[icol,icol]
        a[icol,icol] = 1.0
        for l=1:n
            a[icol,l] = a[icol,l]*piv
        end
        
        for ll=1:n
            if (ll != icol)
                dum = a[ll,icol]
                a[ll,icol] = 0.0
                for l=1:n
                    a[ll,l] = a[ll,l] - a[icol,l]*dum
                end
            end
        end
    end
    
    ## Unscramble Matrix
    
    for l=n:-1:1
        if (indr[l] != indc[l])
            for k = 1:n
                dum = a[k,indr[l]]
                a[k,indr[l]] = a[k,indc[l]]
                a[k,indc[l]] = dum
            end
        end
    end

end 

function vandeven_modal(kk,ngl,p)
    
    ## Constants - ERF
    pe=0.3275911
    a1=0.254829592
    a2=-0.284496736
    a3=1.421413741
    a4=-1.453152027
    a5=1.061405429

    ## Constants - Vandeven
    n=ngl-1
    k=kk-1
    i=2*n/3
    eps=1.0e-10

    if (k <= i)
        x = 0
        return 1
    elseif (k > i && k < n)
        x = Float64(k-i)/Float64(n - i)
        omega = abs(x) - 0.5
        xlog = log(1.0-4.0*omega^2)
        c = 4.0 * omega^2
        diff = abs(x-0.5)
        if (diff < eps)
            square_root = 1
        else
            square_root = sqrt(-xlog/c)
        end
        
        z = 2.0 * sqrt(p) * omega * square_root
        zc = abs(z)
        
        ## ERF
        t = 1.0/(1.0 + pe * zc)
        c = 1.0 - (a1 * t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5) * exp(-zc*zc)
        if (zc < eps)
            c = 0.0
        else
            c = c*z/zc
        end
        return 0.5 * (1.0 - c)
        
    elseif(k == n)
        x=1
        return 0.0
    else
        println(" # problem in Vandeven_modal")
        exit
    end
end
