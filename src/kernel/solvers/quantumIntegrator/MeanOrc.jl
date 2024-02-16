function MeanOrc(Gij, t, N)
    #MEANORC is an oracle function for mean value of g_ij.
    #   MeanOrc is an oracle function that evaluates the mean value of 
    #    G_ij over N knot times. Because of the exponential 
    #    overhead in a Schrodinger simulation of the Quantum Amplitude 
    #    Estimation Algorithm [QAEA], this oracle substitutes for this 
    #    Schrodinger simulation by evaluating the mean value aTrue which is() 
    #    used in QAmpEst as peak in the QAEA probability distribution.
    #
    #   INPUTS: Gij = 1 x N array with values of function g_ij to be averaged
    #                   at N knot times
    #           t = 1 x N array with the N knot times
    #           N = number of knot time in subsubinterval j
    #
    #   OUTPUT:
    #           aTrue = mean value of Gij
    
     temp = 0;      # used to accumulate mean value
     
        #  accumulate mean value
     
        for j = 1: Int(N)
            temp = temp + Gij[j]
        end
     
         aTrue = temp/N;    # RHS is mean value
     
    return aTrue
    end
