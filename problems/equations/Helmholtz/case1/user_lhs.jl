function user_lhs(mesh::St_mesh)

    M = zeros(mesh.npoin,mesh.npoin)
    for ip = 1:mesh.npoin
        M[ip,ip] = 10
    end
    return M
end
