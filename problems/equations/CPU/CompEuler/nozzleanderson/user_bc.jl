function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, qe::SubArray{Float64},::TOTAL)
    
    U1in = 0.0
    U2in = 0.0
    U3in = 0.0
    
    U1out = 0.0
    U2out = 0.0
    U3out = 0.0
    
    ip2 = 2 #this is the 2nd point of the linear grid
    ip3 = 3 #this is the 3rd point of the linear grid
    ipN = npoin_linear #last geometric point of the 1D mesh. The H-O node count starts from this one in the first element.

    xin = 0.0
    Ain = 1.0 + 2.2*(xin - 1.5)^2
    xout = 3.0
    Aout = 1.0 + 2.2*(xout - 1.5)^2
    
    Tin = 1.0
    ρin = 1.0
    pin = ρin*Tin
    γ = 1.4
    γm1 = 0.4
    
    if (tag == "left")
        U1in = Ain
        
        #U2in = 2*q[ip2,2] - q[ip3,2]
        U2in = 2*q[ipN+1,2] - q[ipN+2,2] #0.585
        uin  = U2in/U1in
        U3in = U1in*(Tin/γm1 + 0.5*γ*uin*uin)
        
        qbdy[1] = U1in
        #qbdy[2] = U2in
        qbdy[3] = U3in

       # @info "U2 inner points: " U1in U2in U3in
        
    end

    if (tag == "right")
        pout = 0.6784
            
        #U1out = 2*q[ipN,1] - q[ipN-1,1]
        #U2out = 2*q[ipN,2] - q[ipN-1,2]
        #U3out = 2*q[ipN,3] - q[ipN-1,3]

        #uout = U2out/U1out
        #U3out = pout*Aout/γm1 + 0.5*γ*U2out*uout
        
        #qbdy[1] = U1out
        #qbdy[2] = U2out
        #qbdy[3] = U3out
    end
    
end

function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe::SubArray{Float64},::PERT)
#    if (tag == "free_slip")
      
      if (y<=14950) #(abs(x) < 119500.0 && y<= 14950.0)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
        qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
        qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
      end
      #else 
       # qbdy[2] = 0.0
        #qbdy[3] = 0.0
      #end
      #if (abs(x) > 119500.0 && y < 0.1)
      #  qbdy[2] = 0.0
      #  qbdy[3] = 0.0
      #end
     #@info x,y,nx,ny,qbdy[2],qbdy[3] 
  # return qbdy #, flags
    
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end
