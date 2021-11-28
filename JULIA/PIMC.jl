# Contains all the code that actually performs the PIMC, but not the requisite
# helper functions that they call - obviously not meant to be called on its own

using Printf

# This should probably go in the "action.jl" file
function InitializeDeterminants(Param::Params, Path::Paths)
    for tSlice = 1:Param.nTsl
        for ptcl = 1:Param.nPar
            Path.determinants[tSlice,ptcl] = Determinant(Param, Path, tSlice)
        end
    end
end

function CenterOfMassMove(Param::Params, Path::Paths, ptcl::Int64)
    # Attempts a CoM update, displacing the entire particle worldline
    delta = Param.delta
    shift = Shift(delta) 
    oldAction = 0.0
    newAction = 0.0
    tau = Param.tau
    tauO2 = tau/2

    # Store the postiions on the worldline
    for tSlice = 1:Param.nTsl
        oldAction += ComputeAction(Param, Path,tSlice)
    end
    oldAction *= tauO2

    # Save old potentials to avoid recalculating them
    oldPotentials = copy(Path.potentials)
    if Param.nPar == 1
        # For 1D, determinants don't need updating - translational symmetry for line shift
    else
        # For >1D, definitely need to recompute determinants
    end

    for tSlice = 1:Param.nTsl
        Path.beads[tSlice,ptcl] = Path.beads[tSlice,ptcl] + shift
        UpdatePotential(Path, tSlice, ptcl, Param.lam)
    end

    for tSlice = 1:Param.nTsl
        newAction += ComputeAction(Param, Path,tSlice)
    end
    newAction *= tauO2

    if rand(MersenneTwister()) < exp(-(newAction - oldAction))
        Path.numAcceptCOM += 1
    else
        for tSlice = 1:Param.nTsl
            Path.beads[tSlice,ptcl] = Path.beads[tSlice,ptcl] - shift
        end
        
        Path.potentials = copy(oldPotentials)   # Restore old potentials
    end
end

function StagingMove(Param::Params, Path::Paths, ptcl::Int64)
    #=
    Attempts a staging move, which exactly samples the free-particle propagator
    between two positions.

    See: http://link.aps.org/doi/10.1103/PhysRevB.31.4234

    Note: does not work for periodic boundary conditions.
    =#
    oldAction       = 0.0
    newAction       = 0.0
    edgeExclusion   = 2   # Totally arbitrary choice
        # The length of the stage - must be less than numTimeSlices
    #m               = Param.nTsl - 2*edgeExclusion
    m               = 16
    tau             = Param.tau
    tauO2           = tau/2
    oldBeads        = zeros(m-1)
    oldPotentials   = zeros(m-1)
    oldDeterminant  = zeros(m-1)

    # Choose the start and end of the stage
    alpha_start     = rand(1:Param.nTsl)    # This needs to be inclusive - it is here
    alpha_end       = ModTslice((alpha_start + m), Param.nTsl)

    for a = 1:m-1
        tSlice              = ModTslice(alpha_start + a, Param.nTsl)
        oldBeads[a]         = Path.beads[tSlice,ptcl]
        oldPotentials[a]    = Path.potentials[tSlice,ptcl]
        oldDeterminant[a]   = Path.determinants[tSlice,ptcl]
        oldAction           += tauO2 * ComputeAction(Param, Path, tSlice)
    end

    for a = 1:m-1
        tSlice                      = ModTslice((alpha_start + a), Param.nTsl)
        temp                        = ModTslice((tSlice - 1), Param.nTsl)
        temp == 0 ? tSlicem1 = Param.nTsl : tSlicem1 = temp
        tau1                        = (m - a) * tau
        avex                = (tau1 * Path.beads[tSlicem1,ptcl] + 
                                tau * Path.beads[alpha_end,ptcl]) / (tau + tau1)
        sigma2                      = 2.0 * Param.lam / (1.0 / tau + 1.0 / tau1)
        Path.beads[tSlice,ptcl]     = avex + sqrt(sigma2) * randn()
        UpdatePotential(Path,tSlice, ptcl, Param.lam)
        UpdateDeterminant(Param, Path, tSlice, ptcl)
        newAction                   += tauO2 * ComputeAction(Param, Path,tSlice)
            # See Ceperley about this (Potential) Action, how it relates to the Primitive approximation, extra factors of tau in the "accuracy", etc. In the first few sections.
    end
    
    # Perform the Metropolis step, if we reject, revert the worldline
    if (rand() < exp(-(newAction - oldAction)))
        Path.numAcceptStaging += 1
    else
        for a = 1:m-1
            tSlice = ModTslice((alpha_start + a), Param.nTsl)
            Path.beads[tSlice,ptcl] = oldBeads[a]
            Path.potentials[tSlice,ptcl] = oldPotentials[a]
            Path.determinants[tSlice,ptcl] = oldDeterminant[a]
        end
    end
end

#function Bin(Param::Params, Path::Paths, binArray::Vector{Float64})
function Bin(Param::Params, Path::Paths)
    binLoc = -1

    binArray = zeros(Float64, Param.numMCbins)
    for tSlice = 1:Param.nTsl
        # DON'T FORGET THE "JULIA OFFSET" - +1 BECAUSE JULIA
        binLoc = 1 + trunc(Int, (Path.beads[tSlice,1] - Param.x_a)/Param.binWidth)
        if binLoc < 1
            binLoc = 1
        elseif binLoc > Param.numMCbins
            binLoc = Param.numMCbins
        end
        binArray[binLoc] += 1
    end
    
    return binArray
end

function BinCver(Param::Params, Path::Paths, binArrCount::Vector{Float64})
    binLoc = -1
    for i = 1:Param.numMCbins
        binArrCount[i] = 0
    end

    for tSlice = 1:Param.nTsl
        # DON'T FORGET THE "JULIA OFFSET" - +1 BECAUSE JULIA
        binLoc = 1 + trunc(Int, (Path.beads[tSlice,1] - Param.x_a)/Param.binWidth)
        if binLoc < 1
            binLoc = 1
        elseif binLoc > Param.numMCbins
            binLoc = Param.numMCbins
        end
        binArrCount[binLoc] += 1
    end
end

function UpdateMC(Param::Params, Path::Paths)
    for time = 1:Param.nTsl        
        # pseudo: 
        #   switch: pick a random value
        #       case: CoM move
        #       case: Staging move
        #   end
        for ptcl = 1:rand(1:Param.nPar)
             CenterOfMassMove(Param, Path, ptcl)
        end
        
        for ptcl = 1:rand(1:Param.nPar)
             StagingMove(Param, Path, ptcl)
        end
    end
end

function PIMC(Param::Params, Path::Paths, numSteps::Int64, set::Dict{String, Any})
#= TODO: WHILE I'M AT IT, RE/NAME VARIABLES TO SOMETHING INTELLIGIBLE
# TODO: IMPLEMENT LOGGING OF RESULTS A LA (g)ce-log-### FILE FROM PRODUCTION CODE
# Monte Carlo outline:
#
# for step = 1:numMCsteps
#   for site = 1:numSites
#   TODO: IMPLEMENT THIS RANDOM-SELECTION UPDATE OR LEAVE IT AS-IS?
#       (pick random "site" to perform randomly-chosen update (either CoM and Staging - pick one), do it, iterate)
#       (for the case of this wordline method, pick random line for CoM, random string of beads for Staging)
#           (Check to make sure that this isn't already handled in the CoM and Staging functions (I think it is))
#   end
#   
#   (check - how many steps before measurement?)
#   if (step % observSkip) && (step > equilibriation)
#       (This is to help mitigate the effects of autocorrelation as well as to prevent expensive functions from being computed more than they really need to be)
#       AppropriateAccumulators += ComputeEstimators(arguments)
#       count += 1 (!!! When initializing, don't forget the Julia offset of +1!)
#
#   TODO: IMPLEMENT THIS MC BINNING
#       (Check if the number of accumulated MC steps is the desired amount for result binning - if so, average and then write to disk)
#       if count == numBins
#          (Take the average of the values in the accumulators, write them to disk, then zero the accumulators and the count (technically, `count = 1` is the code) and do it again)
#          numBinsWritten += 1
#       end
#   end
# end
=#
    x1_ave      = 0.0
    x2_ave      = 0.0
    equilSkip   = Param.numEquilibSteps
#    leng        = 1 + trunc(Int,(Param.numSamples - 1) *
#                                Param.numEquilibSteps/Param.observableSkip)
    width       = Param.numMCbins
    energy      = 0.0
    energy2     = 0.0
    ke          = 0.0
    pe          = 0.0
    tau         = Param.tau
    
    distributionArray   = zeros(Float64, width)
    distArrayCount      = zeros(Float64, width)
    stepSize            = (Param.x_b - Param.x_a)/width
        # These are the x-axis bins into with the distribution will be binned
    distrbtnBins                = collect(range(Param.x_a, Param.x_b, step=stepSize))
    deleteat!(distrbtnBins, length(distrbtnBins))   # Stupid range function - it doesn't act 
                                    # like it does in Python. In Julia, it
                                    # includes the endpoint. Thus, remove last.
    if length(distributionArray) != length(distrbtnBins)
        println("Failure! distributionArray = $(length(distributionArray)) while bin = $(length(distrbtnBins))")
        exit()
    end

    # Thanks I hate it
    estDatName  = "ce-estimators-" * Param.baseName
    binDatName  = "ce-lineardensity-" * Param.baseName
    
    WriteHeader(binDatName, "# PIMCID: $(Param.uid)\t#\t#\t#\t#\t#")
    WriteHeader(binDatName, "# ESTINF: dz = $((Param.x_b - Param.x_a)/width) NGRIDSEP = $width\t#\t#\t#")
    
    open(binDatName, "a") do file
        print(file, "#  ")
        for i = 1:length(distrbtnBins)
            @printf(file, "%e   ", distrbtnBins[i])
        end
        print(file, "\n")
    end

    # Initialize the determinants and potentials arrays
    InitializeDeterminants(Param, Path)
    InstantiatePotentials(Param, Path)
    
    # MC iterations
    sampleCount = 0 # Counts number of samples taken
    steps = 1       # Update counter

    # Warmup
    println("Equilibriating the simulation...")
    for steps = ProgressBar(1:equilSkip)
        UpdateMC(Param, Path)
    end

    steps = equilSkip + 1   # Because Julia doesn't "reuse" variables like C/C++ would
    while sampleCount < Param.numSamples
        # The toy code attempts both updates - should I do that here as well?
        UpdateMC(Param, Path)

        if (steps % Param.observableSkip == 0) #&& (steps > equilSkip)  # Should already be handled from above
            sampleCount += 1
            println("Writing binning # $sampleCount / $(Param.numSamples)")
            BinCver(Param, Path, distArrayCount)
            distributionArray   += distArrayCount
            x1_ave              = 0.0
            x2_ave              = 0.0
            energy              = Energy(Param, Path)
            energy2             = energy*energy
            ke                  = Path.KE
            pe                  = Path.PE

            for i = 1:Param.nTsl
                for j = 1:Param.nPar
                    x1_ave += Path.beads[i,j]
                    x2_ave += Path.beads[i,j] * Path.beads[i,j]
                end
            end
            

            WriteHeader(estDatName, "$energy\t$energy2\t$ke\t$pe\t$(x1_ave/Param.nTsl)\t$(x2_ave/Param.nTsl)")
            open(binDatName, "a") do file
                for i = 1:width
                    @printf(file, "  %.8e", distributionArray[i]/(Param.nTsl))
                end
                print(file, "\n")
            end

            for i = 1:width
                distributionArray[i] = 0.0
            end
        end

        steps += 1
    end
end
