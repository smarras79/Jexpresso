include("randQAEA.jl");

using Statistics


function QAmpEst(M, delta, omega)
#QAMPEST Estimates unknown quantum amplitude using QAEA.
#   Function uses Quantum Amplitude Estimation algorithm [QAEA]
#   to estimate the unknown quantum amplitude a = (sin(pi*omega))^2.
#   (See Brassard et al., quant-ph/0005055).
#
#   (i) Function estimates number of runs [TotRuns] needed to insure
#       error probability for final estimate is less than delta.
#       TotRuns calculated in two steps: (a) TempTot = ceil(-8*log(delta))
#       (b) TotRuns equal smallest odd integer greater than TempTot which()
#       insures median calculation returns an element of the Estimates
#       array.
#   (ii) Function generates TotRun estimates for y [=M*omega/pi] &
#        stores them in Estimates[TotRuns] array. Median of this array
#        is estimate of y; amplitude estimate aEstimate = (sin(pi*y/M))^2.
#   (iii) Error on final estimate is O[1/M].
#
#   ALSO: randQAEA.

# true value of unknown amplitude

a = (sin(pi*omega))^2

trueValue = a

# calculate upper bound on estimate error()

UpprBnd = (2*pi/M)*sqrt(a*(1-a)) + (pi/M)^2

# calculate total number of amplitude estimates needed [TotRuns]


TempTot = ceil(-8*log(delta))

if (TempTot%2 .== 0)
TotRuns = TempTot + 1
else
TotRuns = TempTot
end

Estimates = zeros(1,Int(TotRuns))

# start loop to carry out TotRuns simulation runs

for runs = 1:TotRuns
    ###### qaea returns amplitude estimation using algorithm described in review ######
randev = randQAEA(M,omega)   # randQAEA generates random deviate 
           # with probability distribution
           # produced by QAEA
runs = Int(runs)
Estimates[runs] = randev
end

EstimateMedian = median(Estimates)

aEstimate = (sin(pi*EstimateMedian/M))^(2) #step 6 of SI-6C in paper: y from equation is estimateMedian in code (measurement outcome)

error = abs(aEstimate - trueValue)

# no longer interested in this check
#
#if error <= UpprBnd
#message = "Estimate error upper bound satisfied!"
#SuccessFlag = 1
#elseif error > UpprBnd
#message = "Estimate error upper bound violated! "
#SuccessFlag = 0
#end

return [aEstimate, trueValue, error,UpprBnd]

end
