function bathymetry(x,y)

    return -2.0

end

#=function bathymetry(x)

    return -2.0

end=#

#1D parabolic bump Swashes case 1
function bathymetry(x)

    if (x >= 8 && x <= 12)
        Hb = 0.2 - 0.05 * (x - 10)^2
    else
        Hb = 0.0
    end 
    return Hb
end
