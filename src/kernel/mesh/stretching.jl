function stretch_mesh!(mesh,inputs)

    stretch_factor = inputs[:stretch_factor]
    ztop = maximum(mesh.y)
    zsurf = zeros(Float64,mesh.npoin)
    
    for ip = 1:mesh.npoin
        stretching_factor = inputs[:stretch_factor]
        sigma = mesh.y[ip]

        # Normalize the sigma coordinate to the range [0, 1]
        sigma_normalized = sigma / ztop

        # Apply the non-linear power-law stretching function.
        # This maps the uniform [0, 1] grid to a stretched [0, 1] grid where points
        # are concentrated towards 0.
        z_normalized = sigma_normalized ^ stretching_factor

        # Map the stretched, normalized coordinate back to the physical z-domain,
        # which spans from the bottom surface zsurf[ip] to the domain top ztop.
        z = zsurf[ip] + z_normalized * (ztop - zsurf[ip])

        # Update the grid point's vertical position with the new stretched value.
        mesh.y[ip] = z
    end
    
end

function stretch_mesh_3D!(mesh,inputs)  
    
    stretch_factor = inputs[:stretch_factor]
    ztop = maximum(mesh.z)
    zsurf = zeros(Float64,mesh.npoin)

    for ip = 1:mesh.npoin
        stretching_factor = inputs[:stretch_factor]
        sigma = mesh.z[ip]

        # Normalize the sigma coordinate to the range [0, 1]
        sigma_normalized = sigma / ztop

        # Apply the non-linear power-law stretching function.
        # This maps the uniform [0, 1] grid to a stretched [0, 1] grid where points
        # are concentrated towards 0.
        z_normalized = sigma_normalized ^ stretching_factor

        # Map the stretched, normalized coordinate back to the physical z-domain,
        # which spans from the bottom surface zsurf[ip] to the domain top ztop.
        z = zsurf[ip] + z_normalized * (ztop - zsurf[ip])

        # Update the grid point's vertical position with the new stretched value.
        mesh.z[ip] = z
    end
    
end
