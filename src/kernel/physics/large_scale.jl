function read_large_scale!(backend, flist, LST, mesh)
    
    data_out = KernelAbstractions.zeros(backend,TFloat,size(flist,1), Int64(mesh.npoin))
    
    for i=1:size(flist,1)
        data          = read_sounding(flist[1])
        data_reordered = zeros(TFloat, size(data))
        data_reordered[:,1] .= data[:,2]*1000
        data_reordered[:,2] .= data[:,1]
        data_out[i,:] = interpolate_sounding(backend,mesh.npoin,mesh.z,data)
    end

    LST.Rad_cool .= data_out[2,:]
    LST.T_adv    .= data_out[1,:]
    LST.q_adv    .= data_out[3,:]

end

function large_scale_source!(q, qe, S, Rad_cool, T_adv, q_adv,::PERT)
    PhysConst = PhysicalConst{Float64}()
    ρ = q[1] + qe[1]
    S[5] += ρ * (T_adv + Rad_cool) * PhysConst.cp/86400.0
    S[6] += ρ * (q_adv/86400.0)/1000
end
