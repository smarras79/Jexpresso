
function randQAEA(M, omega)
    Momega = M * omega

    SubIntCounter = Int(1)  # SubIntCounter identifies which subinterval is
    #  currently being processed in while-loop below.
    #  Subintervals are numbered 1 -> M, and subinterval
    #  j corresponds to QAEA probability p(j - 1).

    Probs = zeros(1, Int(M)) # stores QAEA prob. dist. p(j-1), with 1<=j<=M

    A = zeros(1, Int(M))     # stores partial sums of p(j-1) from 1 to M terms

    y = range(0, Int(M) - 1, Int(M)) # stores integer values at which p(j-1) is
    # evaluated

    x = y .- Momega     # shifted values appearing in p(y)

    # calculate QAEA probabilities p(j-1), store in Probs(j) with 1<=j<=M
    # note that Probs(1) = p(0), probs(j) = p(j-1), and Probs(M) = p(M-1)

    tempProb = map(sin, pi .* x) ./ map(sin, (pi / M) .* x)

    Probs = (1 / M)^(2) .* (tempProb .* tempProb)

    # accumulate partial sums of p(j-1), 1<=j<=M, store in A, where
    #   A(j) = p(0) + ... + p(j-1)

    A[1] = Probs[1]    # recall Probs(1) = p(0)

    for j = 2:M
        j = Int(j)
        A[Int(j)] = A[Int(j)-1] + Probs[Int(j)]    # recall Probs(j) = p(j-1) 
    end

    # sample uniform deviate u

    u = rand(1)
    #@info u
    # determine which subinterval contains u
    #  1. Loop through subintervals, with SubIntCounter tracking which
    #      subinterval is currently being processed
    #  2. Exit loop when u > A(SubIntCounter) first occurs. Subinterval
    #      labeled by SubIntCounter at exit contains u

    #@info A[SubIntCounter]
    while A[Int(SubIntCounter)] < u[1]
        SubIntCounter = SubIntCounter + 1  # update next subinterval to be
    end                                     #  processed

    # when execution reaches here, Subinterval SubIntCounter contains u.
    #  randev to be returned is SubIntCounter - 1
    randev = SubIntCounter - 1
    #@info randev
    return randev
end