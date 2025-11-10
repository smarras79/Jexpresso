#-----------------------------------------
# 2D Source term function
#-----------------------------------------
function user_source!(S,
                      q, 
                      qe,
                      gradients,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=8,
                      vargs...)

    """
    2D Non-conservative source term for GLM-MHD system:
    Y = (∇·B)(0, B, v·B, v, 0) + (0, 0, ψ(v·∇ψ), 0, v·∇ψ)
    
    Note: In 2D, ∇·B = ∂Bx/∂x + ∂By/∂y (no ∂Bz/∂z term)
    """
    
    # Extract state
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    Bx = q[5]
    By = q[6]
    Bz = q[7]
    ψ  = q[8]
    
    # Primitive velocities
    u = ρu/ρ
    v = ρv/ρ
    w = 0.0  # No z-velocity in pure 2D
    
    # Extract gradients (only x and y derivatives in 2D)
    ∂Bx_∂x = gradients.dBx_dx
    ∂By_∂y = gradients.dBy_dy
    ∂ψ_∂x  = gradients.dψ_dx
    ∂ψ_∂y  = gradients.dψ_dy
    
    # Compute divergence of B (2D: only x and y components)
    divB = ∂Bx_∂x + ∂By_∂y
    
    # Compute v·∇ψ (2D)
    v_dot_gradψ = u*∂ψ_∂x + v*∂ψ_∂y
    
    # Compute v·B (includes out-of-plane Bz!)
    v_dot_B = u*Bx + v*By + w*Bz
    
    # Assemble source term
    S[1] = 0.0
    S[2] = divB * Bx        # x-momentum
    S[3] = divB * By        # y-momentum
    S[4] = divB * v_dot_B + ψ * v_dot_gradψ  # energy
    S[5] = divB * u         # Bx
    S[6] = divB * v         # By
    S[7] = divB * w         # Bz (w=0 in pure 2D, so this is zero)
    S[8] = v_dot_gradψ      # ψ
    
    # Optional: Add damping for ψ
    if hasfield(typeof(PhysConst), :c_p)
        c_h = PhysConst.c_h
        c_p = PhysConst.c_p
        S[8] += -(c_h^2 / c_p^2) * ψ
    end
    
    return nothing
end

function user_source!(S,
                      q, 
                      qe,
                      gradients,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=9
                      vargs...)
    
    """
    Non-conservative source term for GLM-MHD system:
    Y = (∇·B)(0, B, v·B, v, 0) + (0, 0, ψ(v·∇ψ), 0, v·∇ψ)
    
    Note: This must be added to your RHS calculation separately from fluxes!
    """
    PhysConst = PhysicalConst{Float64}()
    
    # Extract state
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    Bx = q[6]
    By = q[7]
    Bz = q[8]
    ψ  = q[9]
    
    # Primitive velocities
    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ
    
    # Extract gradients (these need to be computed in your code!)
    ∂Bx_∂x, ∂By_∂y, ∂Bz_∂z = gradients.dBx_dx, gradients.dBy_dy, gradients.dBz_dz
    ∂ψ_∂x,  ∂ψ_∂y,  ∂ψ_∂z  = gradients.dψ_dx,  gradients.dψ_dy,  gradients.dψ_dz
    
    # Compute divergence of B
    divB = ∂Bx_∂x + ∂By_∂y + ∂Bz_∂z
    
    # Compute v·∇ψ
    v_dot_gradψ = u*∂ψ_∂x + v*∂ψ_∂y + w*∂ψ_∂z
    
    # Compute v·B
    v_dot_B = u*Bx + v*By + w*Bz
    
    # Assemble source term
    S[1] = 0.0
    S[2] = divB * Bx
    S[3] = divB * By
    S[4] = divB * Bz
    S[5] = divB * v_dot_B + ψ * v_dot_gradψ
    S[6] = divB * u
    S[7] = divB * v
    S[8] = divB * w
    S[9] = v_dot_gradψ
    
    # Optional: Add damping for ψ (Dedner et al. 2002)
    # Prevents ψ from growing unbounded
    if hasfield(typeof(PhysConst), :c_p)
        c_h = PhysConst.c_h
        c_p = PhysConst.c_p
        S[9] += -(c_h^2 / c_p^2) * ψ
    end
    
    return nothing
end


function user_source!(S,
                      q, 
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0
   
end

function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ = q[1] #- qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
   
end

function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::NCL,
                      ::AbstractPert; #for NCL() there is no differece between PERT() and TOTAL() in the source
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)
    
    

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -g
    #    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -PhysConst.g
    S[4] = 0.0
    
end

function user_source_gpu(q,qe,x,y,PhysConst, xmax, xmin, ymax, ymin,lpert)

    T = eltype(q)
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    return T(0.0), T(0.0), T(-ρ*PhysConst.g), T(0.0)
end
