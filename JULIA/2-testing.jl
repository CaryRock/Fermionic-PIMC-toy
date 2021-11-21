using Test              # Get the testing going

include("./6-julia_sho.jl")

function TestMain()
# =============================================================================
# Initialization and instantiation
# =============================================================================
    lam = 1 #1/2
    initialT = 0.01
    finalT = 1.00
    numParticles = 1
    numTimeSlices = 20
    numEquilibSteps = 5000
    observableSkip = 500
    numMCSteps = 11*numEquilibSteps + observableSkip
    numProcesses = 10
    temp = finalT*ones(numProcesses)
    tau = collect(range(initialT, finalT, length=numProcesses))
    for index = 1:numProcesses
        tau[index] = 1/(numTimeSlices * finalT)
    end
    binSize = 100
    energyTrace = Float64[]
    potentialTrace = Float64[]
    kineticTrace = Float64[]
    beads = zeros(Float64, numTimeSlices, numParticles)
    for tslice = 1:numTimeSlices
        for ptcl = 1:numParticles
            beads[tslice,ptcl] = 0.5 * (-1.0 + 2.0*rand())
        end
    end
    
    iter = 2
    tslice = 1
    println("tslice = ", tslice)
    println("iter = ", iter)
# =============================================================================
# Testing of complicated functions
# =============================================================================
    testBeads = ones(Float64,numTimeSlices,numParticles)
    mBeads1 = zeros(5,1)
    mBeads2 = zeros(5,1)
    determinants = zeros(Float64,numTimeSlices,numParticles)
    potentials = zeros(Float64,numTimeSlices,numParticles)
    expPots = zeros(Float64,numTimeSlices,numParticles)
    energy = [0 0]
    Path = Paths(numParticles, numTimeSlices, testBeads, lam, tau,
                 determinants, potentials, expPots, 0, 0, 0)
    println("Path.tau[iter=2] = ", Path.tau[iter])

    shift = 0.1
    mBeads1[1,1] =  0.0
    mBeads1[2,1] =  0.5
    mBeads1[3,1] =  0.4
    mBeads1[4,1] = -0.1
    mBeads1[5,1] = -0.5
    for x = 1:5
        mBeads2[x,1] = mBeads1[x,1] + shift
    end

    println("Testing CenterOfMassMove(Path, ptcl, iter): ")
    
    # ComputeAction
    #println("Testing ComputeAction(Path,tslice,iter): ",
    #        @test ComputeAction(Path,tslice,iter) == 

end

if abspath(PROGRAM_FILE) == @__FILE__
    TestMain()
end
