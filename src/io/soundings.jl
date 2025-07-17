function read_sounding(fname)
    data = readdlm(fname,TFloat)
    return data
end

function interpolate_sounding(backend, npoin, z, data)

    n_rows = TInt(size(data,1))
    n_var  = TInt(size(data,2))
    data_out = KernelAbstractions.zeros(backend, TFloat, Int64(npoin),n_var-1)
    if (backend == CPU())
        for ip=1:npoin
            zz = z[ip]
            z1 = minimum(data[:,1])
            i1 = 1
            i2 = n_rows
            z2 = maximum(data[:,1])
            for i=1:n_rows
                z_d = data[i,1]
                if (z_d <= z[ip] && z_d > z1)
                    z1 = z_d
                    i1 = i
                elseif (z_d >= z[ip] && z_d < z2)
                    z2 = z_d
                    i2 = i
                end
            end
            for v=2:n_var
                if (z1 < z2)
                        data_out[ip,v-1] = data[i1,v] + (zz - z1)/(z2-z1)*(data[i2,v] - data[i1,v])
                else
                        data_out[ip,v-1] = data[i1,v]
                end
            end

        end
    else
        data_1 = KernelAbstractions.allocate(backend, TFloat, Int64(n_rows), Int64(n_var))
        KernelAbstractions.copyto!(backend, data_1, data)
        zmin = TFloat(minimum(data[:,1]))
        zmax = TFloat(maximum(data[:,1]))
        k = interpolate_sounding_gpu!(backend)
        k(data_out,z,data_1,n_rows,n_var,zmin,zmax; ndrange = (Int64(npoin)))
    end
    return data_out
end

@kernel function interpolate_sounding_gpu!(data_out,z,data,n_rows,n_var,zmin,zmax)
   ip = @index(Global, Linear)
   z1 = zmin
   z2 = zmax
   T = eltype(data)
   i1 = 1
   i2 = n_rows
   zz = z[ip]
   for i=1:n_rows
        z_d = data[i,1]
        if (z_d <= z[ip] && z_d > z1)
            z1 = z_d
            i1 = i
        elseif (z_d >= z[ip] && z_d < z2)
            z2 = z_d
            i2 = i
        end
    end
    for v=2:n_var
        if (z1 < z2)
            data_out[ip,v-1] = data[i1,v] + T((zz - z1)/(z2-z1))*(data[i2,v] - data[i1,v])
        else
            data_out[ip,v-1] = data[i1,v]
        end
    end
end
