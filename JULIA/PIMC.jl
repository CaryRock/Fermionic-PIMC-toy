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
            # See Ceperley about this (Potential) Action, how it relates to the
            # Primitive approximation, extra factors of tau in the "accuracy", 
            # etc. In the first few sections.
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

function Bin(Param::Params, Path::Paths)
    binLoc = -1

    binArray = zeros(Float64, Param.numSpatialBins)
    for tSlice = 1:Param.nTsl
        # DON'T FORGET THE "JULIA OFFSET" - +1 BECAUSE JULIA
        binLoc = 1 + trunc(Int, (Path.beads[tSlice,1] - Param.x_a)/Param.spatialBinWidth)
        if binLoc < 1
            binLoc = 1
        elseif binLoc > Param.numSpatialBins
            binLoc = Param.numSpatialBins
        end
        binArray[binLoc] += 1
    end
    
    return binArray
end

function SpatialBinCver(Param::Params, Path::Paths, binArrCount::Vector{Float64})
    binLoc = -1
    for i = 1:Param.numSpatialBins
        binArrCount[i] = 0
    end

    for tSlice = 1:Param.nTsl
        # DON'T FORGET THE "JULIA OFFSET" - +1 BECAUSE JULIA
        binLoc = 1 + trunc(Int, (Path.beads[tSlice,1] - Param.x_a)/Param.spatialBinWidth)
        if binLoc < 1
            binLoc = 1
        elseif binLoc > Param.numSpatialBins
            binLoc = Param.numSpatialBins
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
# TODO: IMPLEMENT LOGGING OF RESULTS A LA (g)ce-log-### FILE FROM PRODUCTION CODE
    x1_ave      = 0.0
    x2_ave      = 0.0
    equilSkip   = Param.numEquilibSteps
    width       = Param.numSpatialBins
    energy1     = 0.0
    energy2     = 0.0
    ke          = 0.0
    pe          = 0.0
    tau         = Param.tau
    binSize     = Param.sweepsToBin # It's dumb and I need to fix it

    distributionArray   = zeros(Float64, width)
    distArrayCount      = zeros(Float64, width)
    stepSize            = (Param.x_b - Param.x_a)/width
        # These are the x-axis spatial bins into which the distribution will be binned
    distrbtnBins        = collect(range(Param.x_a, Param.x_b, step=stepSize))
    deleteat!(distrbtnBins, length(distrbtnBins))   # Stupid range function - 
                                    # it doesn't act like it does in Python. 
                                    # In Julia, it includes the endpoint.
                                    # Thus, remove that last point.
    if length(distrbtnBins) != width
        println("Failure! distributionArray = $width while bin = $(length(distrbtnBins))")
        exit()
    end

    estDatName  = "ce-estimator-" * Param.baseName
    if (set["spatial-distribution"])
        binDatName  = "ce-lineardensity-" * Param.baseName
    
        WriteHeader(binDatName, "# PIMCID: $(Param.uid)\t#\t#\t#\t#\t#")
        WriteHeader(binDatName, "# ESTINF: dz = $stepSize NGRIDSEP = $width\t#\t#\t#")
    
        open(binDatName, "a") do file
            print(file, "#  ")
            for i = 1:width
                @printf(file, "%e   ", distrbtnBins[i])
            end
            print(file, "\n")
        end
    end

    # Initialize the determinants and potentials arrays
    InitializeDeterminants(Param, Path)
    InstantiatePotentials(Param, Path)
    
    # MC iterations
    sampleCount = 0 # Counts number of samples taken
    binCount    = 0 # Counts the number of sweeps that have been binned over
    steps       = 1 # Update counter

    # Warmup
    println("Equilibriating the simulation...")
    for steps = ProgressBar(1:equilSkip)
        UpdateMC(Param, Path)
    end

    println("Starting the data collection now...")
    steps = equilSkip + 1   # Because Julia doesn't "reuse" variables like C/C++ would
    while sampleCount < Param.numSamples
        # The toy code attempts both updates - should I do that here as well?
        UpdateMC(Param, Path)

        if (steps % Param.observableSkip == 0)
            ### Start accumulating expectation values, etc.
            binCount += 1

            if (set["spatial-distribution"])
                SpatialBinCver(Param, Path, distArrayCount)
                distributionArray += distArrayCount
            end

            energy1     += Energy(Param, Path)
            energy2     += energy1*energy1
            ke          += Path.KE
            pe          += Path.PE

            for i = 1:Param.nTsl
                for j = 1:Param.nPar
                    x1_ave += Path.beads[i,j]
                    x2_ave += Path.beads[i,j] * Path.beads[i,j]
                end
            end

            if (binCount % binSize == 0)
                ### Write the binned data to file
                sampleCount += 1
                if (sampleCount % 256 == 0)
                    println("Writing sample $sampleCount / $(Param.numSamples)")
                end
                
                energy1 /= binSize
                energy2 /= binSize
                ke      /= binSize 
                pe      /= binSize
                x1_ave  /= (binSize * Param.nTsl)
                x2_ave  /= (binSize * Param.nTsl)

                WriteHeader(estDatName, "$energy1\t$energy2\t$ke\t$pe\t$x1_ave\t$x2_ave")

                if (set["spatial-distribution"])
                    open(binDatName, "a") do file
                        for i = 1:width
                            @printf(file, "  %.8e", distributionArray[i]/(binSize * Param.nTsl))
                        end
                        print(file, "\n")
                    end

                    for i = 1:width
                        distributionArray[i]= 0.0
                        distArrayCount[i]   = 0.0
                    end
                end

                # After recording, zero out the accumulators as well as binCount
                x1_ave  = 0.0
                x2_ave  = 0.0
                energy1 = 0.0
                energy2 = 0.0
            end # end of the binning if()
        end # enf of the observableSkip if()
        
        steps += 1
    end # end of the while(samples) loop

    println("UUID of file(s): $(Param.uid)")
end # end of the PIMC()
