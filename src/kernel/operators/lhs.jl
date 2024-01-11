function build_lhs(L::Array, mesh::St_mesh)
  
    M = user_lhs(mesh)
  
    M .+= L

    return M 

end

function build_lhs(mesh::St_mesh)
  
    M = user_lhs(mesh)

    return M
end
