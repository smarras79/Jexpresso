function initialize(SD::NSD_1D, ET::AdvDiff, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    println(" # Initialize fields for AdvDiff ........................")
    
    qinit = Array{TFloat}(undef, mesh.npoin, 1)
    q     = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, TFloat; neqs=1)
    
    σ = Float64(64.0)
    for iel_g = 1:mesh.nelem
        for i=1:mesh.ngl
            
            ip = mesh.connijk[i,iel_g]
            x  = mesh.x[ip]
            
            #q.qn[ip, 1] = exp(-σ*x*x)
            q.qn[ip, 1] = exp(-200.0*(x - 0.5)^2)
            
        end
    end    
    println(" # Initialize fields for AdvDiff ........................ DONE")
    
    return q
end
