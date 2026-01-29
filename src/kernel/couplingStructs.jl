function couplingAlloc(nrank1, nrank2, T;)
    
    couple = zeros(T, nrank1, nrank2)
    
    return couple
end


#------------------------------------------------------------------------------------
# Receive the grid coordinate from Alya (structured and regular only)
#------------------------------------------------------------------------------------
function distribute_and_count!(
    rem_nx,
    rem_min,
    rem_max,
    ndime,
    nranks2,
    in_my_rank,
    a,
    wrank,
    alya2world)
    
    nx, ny, nz = rem_nx
    nxy        = nx * ny
    nmax       = nxy * nz
    rem_dx     = similar(rem_max)
    
    @assert ndime in (2, 3) "ndime is typically 2 or 3"
    @assert length(rem_min) ≥ ndime && length(rem_nx) ≥ ndime

    r     = mod(nmax, nranks2 - 1)
    npoin = nmax ÷ (nranks2 - 1)

    ri = zeros(Int32, ndime)
    x  = zeros(Float64, ndime)
    
    rem_dx[1:ndime] = (rem_max[1:ndime] .- rem_min[1:ndime])./(rem_nx[1:ndime] .- 1)
    
    @inbounds for ipoin in 1:nmax
        i0   = ipoin - 1
        iz   = i0 ÷ nxy
        remz = i0 - iz * nxy
        iy   = remz ÷ nx
        ix   = mod(remz - iy * nx, nx)

        if ndime ≥ 1; ri[1] = ix; end
        if ndime ≥ 2; ri[2] = iy; end
        if ndime ≥ 3; ri[3] = iz; end

        x[1:ndime] = rem_min[1:ndime] .+ Float32.(ri[1:ndime]) .* rem_dx[1:ndime]

        if in_my_rank
            alya_rank = if ipoin ≤ r * (npoin + 1)
                (ipoin - 1) ÷ (npoin + 1) + 1
            else
                r + (ipoin - r * (npoin + 1) - 1) ÷ npoin + 1
            end

            world_rank = alya2world[alya_rank]
            a[wrank, world_rank] += 1
        end
    end
    println("$a")
    return a
end
