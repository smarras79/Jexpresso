
#Manning-Strickler
function friction(SD::NSD_1D, q, n, npoin; law=1)
    Sf = zeros(npoin)
    C = n^2
    for ip=1:npoin
        if (q[ip,2] > 0.0 && q[ip,1] > 0.1)
            Sf[ip] = C * abs(q[ip,2])*q[ip,2]/(q[ip,1]^(10/3)+1e-16)
        end
    end
    return Sf
end
#=
#Darcy-Weisbach

function friction(SD::NSD_1D, q, f, npoin; law=2)
    @info "Darcy"
    Sf = zeros(npoin)
    C = f/(8*9.81)
    for ip=1:npoin
        Sf[ip] = C * abs(q[ip,2])*q[ip,2]/(q[ip,1]^3)
    end
    return Sf
end         

#Chezy

function friction(SD::NSD_1D, q, C, npoin; law=3)
    @info "Chezy"
    Sf = zeros(npoin)
    Cf = 1/C^2
    for ip=1:npoin
        Sf[ip] = C * abs(q[ip,2])*q[ip,2]/(q[ip,1]^3)
    end
    return Sf
end=#
