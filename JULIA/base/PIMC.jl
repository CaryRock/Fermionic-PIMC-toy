# Contains all the code that actually performs the PIMC, but not all of the 
# requisite helper functions that they call - obviously not meant to be called 
# on its own.

@inbounds function CenterOfMassMove!(Manent::Function, Param::Params, 
                            Path::Paths, ptcl::Int64, rng::MersenneTwister)
    # Attempts a CoM update, displacing the entire particle worldline
    delta = Param.delta
    shift = Shift(rng, delta) 
    oldAction = 0.0::Float64
    newAction = 0.0::Float64


    # Store the postiions on the worldline
    for tSlice = 1:Param.nTsl
        oldAction += ComputeAction(Param, Path, tSlice)
    end
    oldAction *= Param.tau/2.0

    # Save old potentials to avoid recalculating them
    oldPotentials = zeros(Float64, Param.nTsl, Param.nPar)
    for tSlice = 1:Param.nTsl
        for ptcl = 1:Param.nPar
            oldPotentials[tSlice,ptcl] = Path.potentials[tSlice,ptcl]   #@inbounds
        end
    end

    if Param.nPar == 1
        # For a single particle, Path.detemrinants doesn't need updating
    else
        # For multiple particles, Path.determinants does need to be updated. The
        # shifts now also affect interbead, interline spacing
        for tSlice = 1:Param.nTsl
            UpdateManent(Manent, Param, Path, tSlice, ptcl)
        end
    end

    for tSlice = 1:Param.nTsl   # @turbo
        Path.beads[tSlice,ptcl] += shift   # @inbounds
        UpdatePotential(Path, tSlice, ptcl, Param.lam)
    end

    for tSlice = 1:Param.nTsl
        newAction += ComputeAction(Param, Path,tSlice)
    end
    newAction *= Param.tau/2.0

    if rand(rng) < exp(-(newAction - oldAction))    # TODO: FIX THIS SO THAT IT IS ONLY DONE WHEN IT IS SUPPOSED TO
        Path.numAcceptCOM += 1
    else
        for tSlice = 1:Param.nTsl # @turbo
            Path.beads[tSlice,ptcl] -= shift   # @inbounds

            for ptcl = 1:Param.nPar
                Path.potentials[tSlice,ptcl] = oldPotentials[tSlice,ptcl]   # Restore old potentials  # @inbounds
            end
        end
    end
end

@inbounds function StagingMove!(Manent::Function, Param::Params, Path::Paths, 
                                ptcl::Int64, rng::MersenneTwister)
    #=
    Attempts a staging move, which exactly samples the free-particle propagator
    between two positions.

    See: http://link.aps.org/doi/10.1103/PhysRevB.31.4234 

    Note: does not work for periodic boundary conditions.
    =#
    oldAction       = 0.0::Float64
    newAction       = 0.0::Float64
    edgeExclusion   = 1::Int64   # Totally arbitrary choice
        # The length of the stage - must be less than numTimeSlices
    #m               = Param.nTsl - 2*edgeExclusion
    m               = 16::Int64    # TODO: Implement copmmandline switch for this setting
    tau             = Param.tau::Float64
    tauO2           = tau/2.0::Float64
    oldBeads        = zeros(m-1)
    oldPotentials   = zeros(m-1)
    oldDeterminant  = zeros(m-1)

    # Choose the start and end of the stage
    alpha_start     = rand(rng, 1:Param.nTsl)    # This needs to be inclusive - it is here
    alpha_end       = ModTslice((alpha_start + m), Param.nTsl)

    for a = 1:m-1
        tSlice              = ModTslice(alpha_start + a, Param.nTsl)
        oldBeads[a]         = Path.beads[tSlice,ptcl]
        oldPotentials[a]    = Path.potentials[tSlice,ptcl]
        oldDeterminant[a]   = Path.determinants[tSlice,ptcl]
        oldAction           += tauO2 * ComputeAction(Param, Path, tSlice)
    end

    for a = 1:m-1
        tSlice                  = ModTslice((alpha_start + a), Param.nTsl)
        temp                    = ModTslice((tSlice - 1), Param.nTsl)
        temp == 0 ? tSlicem1 = Param.nTsl : tSlicem1 = temp
        tau1                    = (m - a) * tau
        avex                    = (tau1 * Path.beads[tSlicem1,ptcl] + 
                                tau * Path.beads[alpha_end,ptcl]) / (tau + tau1)
        sigma2                  = 2.0 * Param.lam / (1.0 / tau + 1.0 / tau1)
        Path.beads[tSlice,ptcl] = avex + sqrt(sigma2) * randn(rng)
        UpdatePotential(Path, tSlice, ptcl, Param.lam)
        if Param.nPar == 1
        else
            UpdateManent(Manent, Param, Path, tSlice, ptcl)
        end
        newAction                   += tauO2 * ComputeAction(Param, Path,tSlice)
            # See Ceperley about this (Potential) Action, how it relates to the
            # Primitive approximation, extra factors of tau in the "accuracy", 
            # etc. In the first few sections.
    end
    
    # Perform the Metropolis step, if we reject, revert the worldline
    if rand(rng) < exp(-(newAction - oldAction))
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

@inbounds function SpatialBinCver(Param::Params, Path::Paths, binArrCount::Vector{Float64})
    binLoc = -1
    for i = 1:Param.numSpatialBins
        binArrCount[i] = 0
    end

    for tSlice = 1:Param.nTsl
        # DON'T FORGET THE "JULIA OFFSET" - +1 BECAUSE JULIA
        binLoc = 1 + trunc(Int, (Path.beads[tSlice,1] - Param.x_min)/Param.spatialBinWidth)
        if binLoc < 1
            binLoc = 1
        elseif binLoc > Param.numSpatialBins
            binLoc = Param.numSpatialBins
        end
        binArrCount[binLoc] = binArrCount[binLoc] + 1
    end
end

function UpdateMC(Manent::Function, Param::Params, Path::Paths, rng::MersenneTwister)
    for time = 1:Param.nTsl        
        for ptcl = 1:rand(rng, 1:Param.nPar)
            CenterOfMassMove!(Manent, Param, Path, ptcl, rng)
        end
        
        for ptcl = 1:rand(rng, 1:Param.nPar)
            StagingMove!(Manent, Param, Path, ptcl, rng)
        end
    end
end

@inbounds function PIMC(Param::Params, Path::Paths, numSteps::Int64, set::Dict{String, Any}, rng::MersenneTwister)
### Set up the required variables, arrays, and determine what kind of particles
# are being simulated. Also, write out some of the various log files.
    x1_ave      = 0.0::Float64
    x2_ave      = 0.0::Float64
    equilSkip   = Param.numEquilibSteps::Int64
    width       = Param.numSpatialBins::Int64
    energy1     = 0.0::Float64
    energy2     = 0.0::Float64
    ke          = 0.0::Float64
    pe          = 0.0::Float64
    tau         = Param.tau::Float64
    binSize     = Param.sweepsToBin::Int64 # It's dumb and I need to fix it

    distributionArray   = zeros(Float64, width)
    distArrayCount      = zeros(Float64, width)
    stepSize            = (Param.x_max - Param.x_min)/(width - 1)::Int64
        # These are the x-axis spatial bins into which the distribution will be binned
    distrbtnBins        = collect(range(Param.x_min, Param.x_max, step=stepSize))

    println("\nX_max = $(Param.x_max)\nX_Min = $(Param.x_min)")
    println("numSpatialBins = $width")
    println("stepSize = $stepSize\n")

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

### Initialize the determinants and potentials arrays #########################
    function Manent() end
    Manent = WhichManent(Manent, Determinant, Permanent, Boltzmannant, set["bosons"], set["boltzmannons"])
    InstantiateManents(Manent, Param, Path)
    InstantiatePotentials(Param, Path)
    
    # MC iterations
    sampleCount = 0::Int64  # Counts number of samples taken
    binCount    = 0::Int64  # Counts the number of sweeps that have been binned over
    steps       = 1::Int64  # Update counter

    # Warmup
    println("Equilibriating the simulation...")
    for steps = ProgressBar(1:equilSkip)
        UpdateMC(Manent, Param, Path, rng)
    end
    
    # Simulation proper
    println("Starting the data collection now...")
    steps = 1   # Because Julia doesn't "reuse" variables like C/C++ would
    while sampleCount < Param.numSamples
        # The toy code attempts both updates - should I do that here as well?
        UpdateMC(Manent, Param, Path, rng)

        if (steps % Param.observableSkip == 0)
            ### Start accumulating expectation values, etc.
            binCount += 1

            if (set["spatial-distribution"])
                SpatialBinCver(Param, Path, distArrayCount)
                distributionArray += distArrayCount
            end
            
            # Records total Energy of system
            energy1     += Energy(Param, Path)
            energy2     += energy1*energy1
            ke          += Path.KE
            pe          += Path.PE

            # SEE HOW THE PRODUCTION CODE HANDLES MULTIPLE PARTICLES - DOES IT
            # REPORT THEM ALL TOGETHER, DOES IT HANDLE THEM INDIVIDUALLY, THEN 
            # ADD THOSE EXPECTATION VALUES TOGETHER, SOMETHING ELSE...?
            
            # TODO: kludged - this only records the position of the first particle
            for i = 1:Param.nTsl
                x1_ave += Path.beads[i,1]
                x2_ave += Path.beads[i,1] * Path.beads[i,1]
            end

            #=
            for i = 1:Param.nTsl
                for j = 1:Param.nPar
                    x1_ave += Path.beads[i,j] / Param.nPar
                    x2_ave += Path.beads[i,j] * Path.beads[i,j] / Param.nPar
                end
            end
            =#
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
                ke      = 0.0
                pe      = 0.0
                energy1 = 0.0
                energy2 = 0.0
            end # end of the binning if()
        end # end of the observableSkip if()
        
        steps += 1
    end # end of the while(samples) loop

    println("UUID of file(s): $(Param.uid)")
end # end of the PIMC()
