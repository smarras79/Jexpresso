using Quadmath

function filter!(u, params, t, uaux, connijk, connijk_lag, Je, Je_lag, SD::NSD_2D,::TOTAL)
  
  u2uaux!(@view(uaux[:,:]), u, params.neqs, params.mesh.npoin)

  ## Subtract background velocity
  #qv = copy(q)
  uaux[:,2:3] .= uaux[:,2:3] .- params.qe[:,2:3]
  ## store Dimension of MxM object

  ## Loop through the elements

  for e=1:params.mesh.nelem
    for j=1:params.mesh.ngl
      for i=1:params.mesh.ngl
        ip = connijk[e,i,j]
        for m =1:params.neqs
          params.q_t[m,i,j] = uaux[ip,m]
        end
      end
    end
  
  ### Construct local derivatives for prognostic variables
  
    for m=1:params.neqs
    
    ##KSI Derivative
      for i=1:params.mesh.ngl
        for j=1:params.mesh.ngl
          params.q_ti[i,j] = 0.0
          for k=1:params.mesh.ngl
            params.q_ti[i,j] += params.fx[i,k] * params.q_t[m,k,j]
          end
        end
      end


    ## ETA Derivative
      ## this is allocating
      #params.fqf[m,:,:] .= params.q_ti * params.fy_t
      for i=1:params.mesh.ngl
        for j=1:params.mesh.ngl
          params.fqf[m,i,j] = 0.0
          for k=1:params.mesh.ngl
            params.fqf[m,i,j] += params.q_ti[i,k] * params.fy_t[k,j]
          end
        end
      end

    ## ETA Derivative
 

    end
    
  ## Do Numerical Integration

    for j=1:params.mesh.ngl
      for i=1:params.mesh.ngl
        ip = params.mesh.connijk[e,i,j]
        for m=1:params.neqs
          params.b[e,i,j,m] += params.fqf[m,i,j] * params.ω[i]*params.ω[j]*Je[e,i,j]
        end
      end
    end
  end
  
  DSS_rhs!(params.B, params.b, connijk, params.mesh.nelem, params.mesh.ngl, params.neqs, SD, params.AD)
  
  for ieq=1:params.neqs
        divide_by_mass_matrix!(@view(params.B[:,ieq]), params.vaux, params.Minv, params.neqs, params.mesh.npoin, params.AD)
  end
  
          uaux .= params.B
          uaux[:,2:3] .+= params.qe[:,2:3]

  uaux2u!(u, @view(uaux[:,:]), params.neqs, params.mesh.npoin)  
end

function filter!(u, params, t, uaux, connijk, Je, SD::NSD_2D,::PERT; connijk_lag=zeros(TFloat,1,1,1), Je_lag=zeros(TFloat,1,1,1))
 
  u2uaux!(@view(uaux[:,:]), u, params.neqs, params.mesh.npoin)

  #fy_t = transpose(params.fy)
  ## Subtract background velocity
  #qv = copy(q)
  #params.uaux[:,2:4] .= params.uaux[:,2:4] .- params.qe[:,2:4]
  ## store Dimension of MxM object

  ## Loop through the elements

  for e=1:params.mesh.nelem
    for j=1:params.mesh.ngl
      for i=1:params.mesh.ngl
        ip = connijk[e,i,j]
        for m =1:params.neqs
          params.q_t[m,i,j] = uaux[ip,m]
        end
      end
    end
 
  ### Construct local derivatives for prognostic variables
   ### this section accouns for 1/3 of the allocations and more than half in terms of storage size 
   ##(159.84 k allocations: 22.544 MiB) current function total, killed 1/3 of allocations thanks to loop unroll
    for m=1:params.neqs
      #this loop unroll works well for both matmuls allocations now: (108.00 k allocations: 9.888 MiB)
      for i=1:params.mesh.ngl 
        for j=1:params.mesh.ngl
          params.q_ti[i,j] = 0.0
          for k=1:params.mesh.ngl
            params.q_ti[i,j] += params.fx[i,k] * params.q_t[m,k,j]
          end
        end
      end


    ## ETA Derivative
      ## this is allocating
      #params.fqf[m,:,:] .= params.q_ti * params.fy_t
      for i=1:params.mesh.ngl
        for j=1:params.mesh.ngl
          params.fqf[m,i,j] = 0.0
          for k=1:params.mesh.ngl
            params.fqf[m,i,j] += params.q_ti[i,k] * params.fy_t[k,j]
          end
        end
      end

    end
   
  ## Do Numerical Integration

    for j=1:params.mesh.ngl
      for i=1:params.mesh.ngl
        for m=1:params.neqs
          params.b[e,i,j,m] += params.fqf[m,i,j] * params.ω[i]*params.ω[j]*Je[e,i,j]
        end
      end
    end
  end

  DSS_rhs!(params.B, params.b, connijk, params.mesh.nelem, params.mesh.ngl, params.neqs, SD, params.AD)

  if (params.laguerre)

    for e=1:params.mesh.nelem_semi_inf
      for j=1:params.mesh.ngr
        for i=1:params.mesh.ngl
          ip = connijk_lag[e,i,j]
          for m =1:params.neqs
            params.q_t_lag[m,i,j] = uaux[ip,m]
          end
        end
      end

  ### Construct local derivatives for prognostic variables
   ### this section accouns for 1/3 of the allocations and more than half in terms of storage size
   ##(159.84 k allocations: 22.544 MiB) current function total, killed 1/3 of allocations thanks to loop unroll
      for m=1:params.neqs
        #this loop unroll works well for both matmuls allocations now: (108.00 k allocations: 9.888 MiB)
        for i=1:params.mesh.ngl
          for j=1:params.mesh.ngr
            params.q_ti_lag[i,j] = 0.0
            for k=1:params.mesh.ngl
              params.q_ti_lag[i,j] += params.fx[i,k] * params.q_t_lag[m,k,j]
            end
          end
        end


    ## ETA Derivative
      ## this is allocating
      #params.fqf[m,:,:] .= params.q_ti * params.fy_t
        for i=1:params.mesh.ngl
          for j=1:params.mesh.ngr
            params.fqf_lag[m,i,j] = 0.0
            for k=1:params.mesh.ngr
              params.fqf_lag[m,i,j] += params.q_ti_lag[i,k] * params.fy_t_lag[k,j]
              #if (k == j)
               # params.fqf_lag[m,i,j] += params.q_ti_lag[i,k] * 1.0
              #else
              #  params.fqf_lag[m,i,j] += params.q_ti_lag[i,k] * 0.0
              #end
            end
          end
        end

      end

      for j=1:params.mesh.ngr
        for i=1:params.mesh.ngl
          for m=1:params.neqs
            params.b_lag[e,i,j,m] += params.fqf_lag[m,i,j] * params.ω[i]*params.ω_lag[j]*Je_lag[e,i,j]
          end
        end
      end
    end

    DSS_rhs_laguerre!(params.B_lag, params.b_lag, connijk_lag, params.mesh.nelem_semi_inf, params.mesh.ngl, params.mesh.ngr, params.neqs, SD, params.AD)
    #for ip=1:params.mesh.npoin
      #if !(ip in params.mesh.poin_in_bdy_edge)
    params.B .+= params.B_lag
        
      #else
        #if (ip in params.mesh.poin_in_bdy_edge && params.mesh.y[ip] > 14000.0 && abs(params.mesh.x[ip]) < 10000.0)
        #@info  t, params.B[ip,:], params.B_lag[ip,:],ip, params.mesh.x[ip],params.mesh.y[ip]
        #end
      #end
    #end
  end

  #@info "before div"
  #@info params.B[3247,:]

  for ieq=1:params.neqs
       divide_by_mass_matrix!(@view(params.B[:,ieq]), params.vaux, params.Minv, params.neqs, params.mesh.npoin, params.AD)
  end
  #@info "after div"
  #@info params.B[3247,:]
  #@info "before filtering"
  #@info params.uaux[3247,:]
  uaux .= params.B

  #=if (params.laguerre)
  
    @time uaux .= params.B

  end=#


  uaux2u!(u, @view(uaux[:,:]), params.neqs, params.mesh.npoin)
end

function filter!(u, params, t, uaux, connijk, Je, SD::NSD_3D,::PERT; connijk_lag=zeros(TFloat,1,1,1,1), Je_lag=zeros(TFloat,1,1,1,1))

  u2uaux!(@view(uaux[:,:]), u, params.neqs, params.mesh.npoin)

  #fy_t = transpose(params.fy)
  ## Subtract background velocity
  #qv = copy(q)
  #params.uaux[:,2:4] .= params.uaux[:,2:4] .- params.qe[:,2:4]
  ## store Dimension of MxM object

  ## Loop through the elements

  for e=1:params.mesh.nelem
    for k=1:params.mesh.ngl
      for j=1:params.mesh.ngl
        for i=1:params.mesh.ngl
          ip = connijk[e,i,j,k]
          for m =1:params.neqs
            params.q_t[m,i,j,k] = uaux[ip,m]
          end
        end
      end
    end

  ### Construct local derivatives for prognostic variables
   ### this section accouns for 1/3 of the allocations and more than half in terms of storage size 
   ##(159.84 k allocations: 22.544 MiB) current function total, killed 1/3 of allocations thanks to loop unroll
    for m=1:params.neqs
      #this loop unroll works well for both matmuls allocations now: (108.00 k allocations: 9.888 MiB)
      for i=1:params.mesh.ngl
        for j=1:params.mesh.ngl
          for k=1:params.mesh.ngl
            params.q_ti[i,j,k] = 0.0
            for l=1:params.mesh.ngl
              params.q_ti[i,j,k] += params.fx[i,l] * params.q_t[m,l,j,k]
            end
          end
        end
      end


    ## ETA Derivative
      ## this is very likely wrong, work out on paper
      #params.fqf[m,:,:] .= params.q_ti * params.fy_t
      for i=1:params.mesh.ngl
        for j=1:params.mesh.ngl
          for k=1:params.mesh.ngl
            params.q_tij[i,j,k] = 0.0
            for l=1:params.mesh.ngl
                params.q_tij[i,j,k] += params.q_ti[i,l,k] * params.fy_t[l,j]
            end
          end
        end
      end

      for i=1:params.mesh.ngl
          for j=1:params.mesh.ngl
              for k =1:params.mesh.ngl
                  params.fqf[m,i,j,k] = 0.0
                  for l=1:params.mesh.ngl
                      params.fqf[m,i,j,k] += params.q_tij[i,j,l] * params.fz_t[l,k]
                  end
              end
          end
      end

    end

  ## Do Numerical Integration

    for j=1:params.mesh.ngl
      for i=1:params.mesh.ngl
          for k=1:params.mesh.ngl
            for m=1:params.neqs
                params.b[e,i,j,k,m] += params.fqf[m,i,j,k] * params.ω[i]*params.ω[j]*params.ω[k]*Je[e,i,j,k]
            end
        end
      end
    end
  end

  DSS_rhs!(params.B, params.b, connijk, params.mesh.nelem, params.mesh.ngl, params.neqs, SD, params.AD)

  if (params.laguerre)

    for e=1:params.mesh.nelem_semi_inf
      for j=1:params.mesh.ngr
        for i=1:params.mesh.ngl
          ip = connijk_lag[e,i,j]
          for m =1:params.neqs
            params.q_t_lag[m,i,j] = uaux[ip,m]
          end
        end
      end

  ### Construct local derivatives for prognostic variables
   ### this section accouns for 1/3 of the allocations and more than half in terms of storage size
   ##(159.84 k allocations: 22.544 MiB) current function total, killed 1/3 of allocations thanks to loop unroll
      for m=1:params.neqs
        #this loop unroll works well for both matmuls allocations now: (108.00 k allocations: 9.888 MiB)
        for i=1:params.mesh.ngl
          for j=1:params.mesh.ngr
            params.q_ti_lag[i,j] = 0.0
            for k=1:params.mesh.ngl
              params.q_ti_lag[i,j] += params.fx[i,k] * params.q_t_lag[m,k,j]
            end
          end
        end

        for i=1:params.mesh.ngl
          for j=1:params.mesh.ngr
            params.fqf_lag[m,i,j] = 0.0
            for k=1:params.mesh.ngr
              params.fqf_lag[m,i,j] += params.q_ti_lag[i,k] * params.fy_t_lag[k,j]
              #if (k == j)
               # params.fqf_lag[m,i,j] += params.q_ti_lag[i,k] * 1.0
              #else
              #  params.fqf_lag[m,i,j] += params.q_ti_lag[i,k] * 0.0
              #end
            end
          end
        end

      end

      for j=1:params.mesh.ngr
        for i=1:params.mesh.ngl
          for m=1:params.neqs
            params.b_lag[e,i,j,m] += params.fqf_lag[m,i,j] * params.ω[i]*params.ω_lag[j]*Je_lag[e,i,j]
          end
        end
      end
    end

    DSS_rhs_laguerre!(params.B_lag, params.b_lag, connijk_lag, params.mesh.nelem_semi_inf, params.mesh.ngl, params.mesh.ngr, params.neqs, SD, params.AD)
    #for ip=1:params.mesh.npoin
      #if !(ip in params.mesh.poin_in_bdy_edge)
    params.B .+= params.B_lag

      #else
        #if (ip in params.mesh.poin_in_bdy_edge && params.mesh.y[ip] > 14000.0 && abs(params.mesh.x[ip]) < 10000.0)
        #@info  t, params.B[ip,:], params.B_lag[ip,:],ip, params.mesh.x[ip],params.mesh.y[ip]
        #end
      #end
    #end
  end

  #@info "before div"
  #@info params.B[3247,:]

  for ieq=1:params.neqs
       divide_by_mass_matrix!(@view(params.B[:,ieq]), params.vaux, params.Minv, params.neqs, params.mesh.npoin, params.AD)
  end

  uaux[:,params.neqs] .= params.B[:,params.neqs]

  #=if (params.laguerre)

    @time uaux .= params.B

  end=#


  uaux2u!(u, @view(uaux[:,:]), params.neqs, params.mesh.npoin)
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

@kernel function filter_gpu_3d!(u, B, fx, fy_t, fz_t, Je, ω_x, ω_y, ω_z, Minv, connijk, n_x, n_y, n_z, neqs)
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
        @inbounds q_t[i,j,k] = u[ip,m]
        @synchronize()
        q_ti[i,j,k] = zero(eltype(u))
        for l=1:n_x
           @inbounds q_ti[i,j,k] += fx[i,l] * q_t[l,j,k]
        end

        @synchronize()
        q_tij[i,j,k] = zero(eltype(u))
        for n=1:n_z
            for l=1:n_y
                @inbounds q_tij[i,j,n] += q_ti[i,l,n] * fy_t[l,j]
            end
        end

        @synchronize()
        fqf[i,j,k] = zero(eltype(u))
        for l=1:n_z
            @inbounds fqf[i,j,k] += q_ti[i,j,k] * fy_z[l,k]
        end


        @inbounds KernelAbstractions.@atomic B[ip,m] += ω_x[i]*ω_y[j]*ω_z*Je[ie,i,j,k]*fqf[i,j,k]
    end 
end


function init_filter(nop,xgl,mu_x,mesh,inputs)

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
    @info "Legendre filter"
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
    @info "Laguerre filter"
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
    @info "Error in GAUJORDF in FILTER INIT"
    @info "ierr", ierr
    exit
  end
  
  ## Compute Boyd-Vandeven (ERF-LOG) Transfer function
  filter_type = inputs[:filter_type]
  if (filter_type == "erf")   
    @info "erf filtering on"
    for k=1:nop+1
      # Boyd filter
      weight[k] = vandeven_modal(k,nop+1,erf_order)
    end
  elseif (filter_type == "quad")
    @info "quadratic filtering on"
    mode_filter = floor(quad_order)   
    k0 = Int64(nop+1 - mode_filter)
    xmode2 = mode_filter*mode_filter
    weight .= 1
    for k=k0+1:nop+1
      amp = quad_alpha*(k-k0)*(k-k0)/(xmode2)
      weight[k] = 1.0 - amp
    end
  elseif (filter_type == "exp")
    @info "exponential filtering on"
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
    @info "problem in Vandeven_modal"
    exit
  end
end
