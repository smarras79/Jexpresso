function Set_XGrid( x_min, x_max, Tot_X_Pts )
    #SET_XGRID creates uniform grid on [x_min, x_max] with Tot_X_Pts points
    #   Set_XGrid creates uniform grid on interval [x_min, x_max] with
    #       Tot_X_Pts points. 
    #
    #   INPUTS:
    #           x_min = initial point of interval
    #           x_max = final point of interval
    #           Tot_X_Pts = total number of grid points created
    #
    #   OUTPUTS:
    #           x = 1 x Tot_X_PTs array with location of grid points
    #           Del_x = width of each grid subinterval 

    # calculate Del_x
    Del_x = (x_max - x_min)/(Tot_X_Pts - 1)

    # initialize array x
    x = zeros(1,Tot_X_Pts)

    for i = 1:Tot_X_Pts
        if i .== 1
            x[i] = x_min
        else()
            x[i] = x[i-1] + Del_x
        end
    end

    return  x, Del_x 
end