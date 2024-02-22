function user_bc_dirichlet!(ip,
                            q::SubArray{Float64},
                            mesh::St_mesh,
                            t::AbstractFloat,
                            tag::String,
                            qbdy::AbstractArray,                            
                            qe::SubArray{Float64},
                            AD, ::TOTAL)

    if AD == FD()
        U1in = 0.0
        U2in = 0.0
        U3in = 0.0
        
        U1out = 0.0
        U2out = 0.0
        U3out = 0.0
        
        ip2 = 2 #this is the 2nd point of the linear grid
        ip3 = 3 #this is the 3rd point of the linear grid
        ipN = length(q[:,1]) #last geometric point of the 1D mesh. The H-O node count starts from this one in the first element.
        
        xin = 0.0
        Ain = 1.0 + 2.2*(xin - 1.5)^2
        xout = 3.0
        Aout = 1.0 + 2.2*(xout - 1.5)^2
        
        Tin = 1.0
        ρin = 1.0
        pin = ρin*Tin
        mass_flow = 0.59
        uin = mass_flow/(ρin*Ain)
        
        γ = 1.4
        γm1 = 0.4
        
        lshock = false #Notice, only try shock if you have some artificial diffusion implemented
        
        if (tag == "left")
            
            U1in = Ain*ρin
            U2in = 2*q[ip2,2] - q[ip3,2]
            uin  = U2in/U1in
            U3in = U1in*(1/γm1 + 0.5*γ*(U2in/U1in)^2)
            qbdy[1] = U1in
            qbdy[2] = U2in
            qbdy[3] = U3in

            # @info "U2 inner points: " U1in U2in U3in
            
        end

        if (tag == "right")
            pout = 0.6784
            
            U1out = 2*q[ipN-1,1] - q[ipN-2,1]
            U2out = 2*q[ipN-1,2] - q[ipN-2,2]
            U3out = 2*q[ipN-1,3] - q[ipN-2,3]
            
            qbdy[1] = U1out
            qbdy[2] = U2out
            qbdy[3] = U3out
        end

    else
      
        #
        # ContGal
        #
        U1in = 0.0
        U2in = 0.0
        U3in = 0.0
        
        U1out = 0.0
        U2out = 0.0
        U3out = 0.0
               
        
        xin = 0.0
        Ain = 1.0 + 2.2*(xin - 1.5)^2
        xout = 3.0
        Aout = 1.0 + 2.2*(xout - 1.5)^2
        
        Tin = 1.0
        ρin = 1.0
        pin = ρin*Tin
        mass_flow = 0.59
        uin = mass_flow/(ρin*Ain)
        
        γ = 1.4
        γm1 = 0.4
        
        lshock = false #Notice, only try shock if you have some artificial diffusion implemented
        
        if (tag == "left")
            
            U1in = Ain
            
            U2in = je_spline_1D_left(mesh, q, mesh.xmin; ielem=1, ivar=2)
            
            uin  = U2in/U1in
            U3in = U1in*(1/γm1 + 0.5*γ*(U2in/U1in)^2)

            qbdy[1] = U1in
            qbdy[2] = U2in
            qbdy[3] = U3in
           
        elseif (tag == "right")
            pout = 0.6784        

            U1out = je_spline_1D_right(mesh, q, mesh.xmax; ielem=mesh.nelem, ivar=1)
            #U1out = 2*q[ipN-1,1] - q[ipN-2,1]
            
            U2out = je_spline_1D_right(mesh, q, mesh.xmax; ielem=mesh.nelem, ivar=2)
            #U2out = 2*q[ipN-1,2] - q[ipN-2,2]
            
            U3out = je_spline_1D_right(mesh, q, mesh.xmax; ielem=mesh.nelem, ivar=3)
            #U3out = 2*q[ipN-1,3] - q[ipN-2,3]
            
            qbdy[1] = U1out
            qbdy[2] = U2out
            qbdy[3] = U3out
        end
        
    end
        
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function je_spline_1D_left(mesh, q, xbdy; ielem=1, ivar=1)
    
    @info "LEFT bdy"
    ip1 = mesh.connijk[ielem, 1]
    ip2 = mesh.connijk[ielem, 2]
    ip3 = mesh.connijk[ielem, 3]
    ip4 = mesh.connijk[ielem, 4]
    ipn = mesh.connijk[ielem, mesh.ngl]

    xloc    = zeros(mesh.ngl-1)
    xloc[1] = mesh.x[ip2]
    xloc[2] = mesh.x[ip3]
    xloc[3] = mesh.x[ip4]
    xloc[4] = mesh.x[ipn]
    
    qloc    = zeros(mesh.ngl-1)
    qloc[1] = q[ip2,ivar]
    qloc[2] = q[ip3,ivar]
    qloc[3] = q[ip4,ivar]
    qloc[4] = q[ipn,ivar]
    
    
    @info xloc[1] xloc[2] xloc[3] xloc[4]
    @info qloc[1] qloc[2] qloc[3] qloc[4]
    
    nodes = (xloc,)

    qi = linear_interpolation(xloc, qloc, extrapolation_bc=Line())
    
    return qi(xbdy)

end

function je_spline_1D_right(mesh, q, xbdy; ielem=mesh.nelem, ivar=1)
    
    @info "RIGHT bdy"
    ip1 = mesh.connijk[ielem, 1]
    ip2 = mesh.connijk[ielem, 2]
    ip3 = mesh.connijk[ielem, 3]
    ip4 = mesh.connijk[ielem, 4]
    ipn = mesh.connijk[ielem, mesh.ngl]

    xloc = zeros(mesh.ngl-1)
    xloc[1] = mesh.x[ip1]
    xloc[2] = mesh.x[ip2]
    xloc[3] = mesh.x[ip3]
    xloc[4] = mesh.x[ip4]
    
    qloc = zeros(mesh.ngl-1)
    qloc[1] = q[ip1,ivar]
    qloc[2] = q[ip2,ivar]
    qloc[3] = q[ip3,ivar]
    qloc[4] = q[ip4,ivar]
    
    @info xloc[1] xloc[2] xloc[3] xloc[4]
    @info qloc[1] qloc[2] qloc[3] qloc[4]
    
    nodes = (xloc,)

    qi = linear_interpolation(xloc, qloc, extrapolation_bc=Line())

    @info xbdy
    @info qi(xbdy)
    return qi(xbdy)

end
