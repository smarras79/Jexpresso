# Rayleigh sponge, protocol 1.3: relax POTENTIAL TEMPERATURE and the three WIND
# COMPONENTS, S = -c*rho*(u - u_e), -c*rho*(theta - theta_e). The old form relaxed
# rho*u and rho*theta with rho free, i.e. pulled theta towards rho_e*theta_e/rho:
# a compressed parcel was cooled and sank harder. With the mountain waves of the
# 1000 m ridge in the sponge that grew a disturbance at 3-4.5 km (w rms 2.06 at
# 3.5 km, rising; alpha 1/60 made it worse) that the primitive form removes
# (0.34, flat; jobs 1401354 / 1402636 / 1403092). Read once at include time.
# DBG_SPONGE_FORM=conserved restores the old form.
const LESICP_SPONGE_PRIMITIVE = get(ENV, "DBG_SPONGE_FORM", "primitive") == "primitive"

function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=0.0,
                      y=0.0,
                      z=0.0,
                      xmin=0.0, xmax=0.0,
                      ymin=0.0, ymax=0.0,
                      zmin=0.0, zmax=0.0)

    PhysConst = PhysicalConst{Float64}()
    
    #--------------
    # S(q(x)) = -ρg
    #--------------
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = -ρ*PhysConst.g
    S[5] = 0.0

    #--------------
    # Coriolis
    #--------------
    lcoriolis = true
    lgeostrophic = true
    
    #--------------
    # SPONGE
    #--------------
        zs = 3000.0
    	xr = 0.0
    	xl = 0.0
    	α  = 1.0/300.0
	if (z >= zs)#nsponge_points * dsy) #&& dbl >= 0.0)
            betay_coe = 0.5*α*(1.0-cospi((z - zs)/(zmax - zs)))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
	else
		betay_coe = 0.0
	end
    	ctop= 1.0*betay_coe

	#if (x >= xr)#nsponge_points * dsy) #&& dbl >= 0.0)
	        #    betaxr_coe =  sinpi(0.5*(x-xr)/(xmax-xr))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    	#else
    	betaxr_coe = 0.0
    	#end
    
       #if (x <= xl)#nsponge_points * dsy) #&& dbl >= 0.0)
       #    betaxl_coe =  sinpi(0.5*(xl-x)/(xl-xmin))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
        #else
         betaxl_coe = 0.0
        #end
    
        cxr = 0.0*betaxr_coe
        cxl = 0.0*betaxl_coe
        cyr = 0.0
        cyl = 0.0
        cs  = 1.0 - (1.0 - ctop)*(1.0 - cxr)*(1.0 - cxl)*(1.0 - cyr)*(1.0 - cyl)

        #@info "β x: " ctop,cxr,cxl,cs, zs, y, x, ymin, ymax, dsy, dbl
        #S[1] -= (cs)*(q[1]-qe[1])
        if LESICP_SPONGE_PRIMITIVE
            ρe = qe[1]
            S[2] -= cs*ρ*(q[2]/ρ - qe[2]/ρe)
            S[3] -= cs*ρ*(q[3]/ρ - qe[3]/ρe)
            S[4] -= cs*ρ*(q[4]/ρ - qe[4]/ρe)
            S[5] -= cs*ρ*(q[5]/ρ - qe[5]/ρe)
        else
            S[2] -= cs*(q[2]-qe[2])
            S[3] -= cs*(q[3]-qe[3])
            S[4] -= cs*(q[4]-qe[4])
            S[5] -= cs*(q[5]-qe[5])
        end

    
    #Coriolis & geostrophic wind    
    if lcoriolis == true
        f = 1.0e-4
        u_vel = q[2]
        v_vel = q[3]
        S[2] += f * v_vel
        S[3] -= f * u_vel

        if lgeostrophic == true
            # Protocol section 1.7: u_geo, v_geo are the FREE-ATMOSPHERE wind of
            # this case (u10_ridge100), a constant -- not the local sounding
            # wind, which is 0 below 1000 m in u10_ridge1000.
            U_geo = 10.0
            V_geo = 0.0
            S[2] -= q[1] * f * V_geo
            S[3] += q[1] * f * U_geo
        end
    end
        
    return  S
end
