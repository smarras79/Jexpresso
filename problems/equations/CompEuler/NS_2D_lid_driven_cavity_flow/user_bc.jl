
function user_bc_dirichlet!(x, y, t, tag, velocity)

    if (tag == "LID")
        velocity[1] = 1
        velocity[2] = 0
    else
        velocity[:] .= 0
    end

end

