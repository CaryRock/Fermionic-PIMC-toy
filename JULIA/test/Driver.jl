# Notes
# * The Max Graves thesis has the derivation of the virial implementation of the
# KE estimator - unnecessary for this code, but will be interesting/useful later

using Dates
using Random
using ArgParse
using Printf
#using LoopVectorization
using ProgressBars  # Entirely vain - just gives the pretty progress bar
using UUIDs         # I sort-of despise this language. I really kind of do.

#using .Threads

using TimerOutputs  # Optimization profiling

struct Params
    nPar::Int64             # Number of particles
    nTsl::Int64             # Number of time slices
    lam::Float64            # System parameters - hbar^2/(2mk_B) = 1/2
    tau::Float64            # beta/J (beta/M in Ceperley-ese)
    x_min::Float64            # Arbitrary left edge
    x_max::Float64            # Arbitrary right edge
    delta::Float64          # Width of possible shift for a bead
    numHistBins::Int64      # Number of bins for histogram
    sweepsToBin::Int64
    # The following two should be specific to relevant estimators, not global
    spatialBinWidth::Float64# Width of bins for binning the distribution
    numSpatialBins::Int64   # Number of bins for MC binning
    numEquilibSteps::Int64  # Number of steps to skip for equilibriation
    observableSkip::Int64   # Number of MC steps to skip between observations
    numSamples::Int64       # Sets # of MC setps in total
    baseName::String        # Base name of data file(s)
    uid::UUID               # UUID
end

mutable struct Paths
    beads::Array{Float64,}
    determinants::Array{Float64}
    potentials::Array{Float64}
    KE::Float64
    PE::Float64
    numAcceptCOM::Float64
    numAcceptStaging::Float64
end

function parse_commandline()
    set = ArgParseSettings()

    @add_arg_table set begin
### Required values ###########################################################
        "--temp", "-T"
            help = "Temperature in K"
            required = true
            arg_type = Float64
            default = 1.0
        "--numPart", "-N"
            help = "Number of particles"
            required = true
            arg_type = Int64
            default = 1
        "--tau", "-t"
            help = "delta tau of particle - 1/(temperature * number of time slices). Defineds the number of time slices."
            required = true
            arg_type = Float64
            default = 0.02      # For T = 1 => J = 50
#        "--numTime", "-J"
#            help = "Number of time slices"
#            required = true
#            arg_type = Int64
        "--numEquil", "-E"
            help = "Number of steps to equilibriate over"
            required = true
            arg_type = Int64
            default = 32768     # 2^15
        "--obsSkip", "-O"
            help = "Number of steps between observations"
            required = true
            arg_type = Int64
            default = 1
        "--numSamples", "-S"
            help = "Number of desired samples to take"
            required = true
            arg_type = Int64
            default = 131072    # 2^17
### Options for Setting Variables #############################################
        "--sweepsToBin", "-B"
            help = "Sets the desired number of MC sweeps to bin over"
            arg_type = Int64
            default = 50
        "--spatialBinWidth", "-w"
            help = "Sets the desired width of bins"
            arg_type = Float64
            default = 0.05
        "--delta", "-d"
            help = "Sets the delta for the CoM move"
            arg_type = Float64
            default = 1.0
        "--Xmin"
            help = "Sets the minimum value for X - can be negative"
            arg_type = Float64
            default = -5.0
        "--Xmax"
            help = "Sets the maximum value for X - can be negative"
            arg_type = Float64
            default = 5.0
        "--bosons", "-b"
            help = "Informs the simulation to use bosons instead of fermions."
            action = :store_true
        "--boltzmannons", "-z"
            help = "Informs the simulation to use \"boltzmannons\" instead of fermions."
            action = :store_true
### Optional Output Flags #####################################################
        "--rawX1"
            help = "Reporti raw <x> - WARNING: GENERATES MUCH DATA"
            action = :store_true
        "--rawX2"
            help = "Report raw <x^2> - WARNING: GENERATES MUCH DATA"
            action = :store_true
        "--binX1"
            help = "Report binned <x> value"
            action = :store_true
        "--binX2"
            help = "Report binned <x^2> value"
            action = :store_true
        "--spatial-distribution"
            help = "Record spatial (1D) distribution to go along with the run."
            action = :store_true
    end
    
    return parse_args(set)
end

# The external potential to use
# If more start being added, this should go into a separate file with the other
# potentials
@inline function ExtPotential(lam::Float64, R::Float64)
    return 1/(4 * lam) * R * R
end

# Computes the modular time slice to use
# +1 because Julia is 1-indexed, not 0-indexed like Python or C
@inline function ModTslice(tSlice::Int64, numSlices::Int64)
    return (( abs(tSlice) - 1 + numSlices ) % numSlices ) + 1
end

# Flags whether or not to truncate the calculation
@inline function CutOff(beadA::Float64, beadB::Float64, cutoff::Float64=10^-4)
    return ( abs(beadA - beadB) < cutoff )
end

# Controls the width of the possible shift in position for a bead
@inline function Shift(rng::MersenneTwister, delta::Float64  = 1.0)
    return delta * (-1.0 + 2.0 * rand(rng))
end

# Writes the header of the given file, filename, with the given string
# Technically, it only writes a single string at a time, so it can be used for 
# writing single lines
function WriteHeader(fileName::String, string)
    open(fileName, "a") do file
        println(file, string)
    end
end

# Writes out a whole matrix work of data (data) to a given file (filename)
function WriteOutMat(fileName::String, data, xLim::Int64 = length(data), yLim::Int64 = length(data[1]))
    open(fileName, "a") do file
        for t = 1:xLim
            for p = 1:yLim
                println(file, data[t,p])
            end
        end
    end
end

# Writes out the first part of the log file - the simulation parameters part
function WriteLogFileParameters(fileName::String, Param::Params, Path::Paths,
        numMCsteps::Int64, set::Dict{String, Any}, invocation::String)
    open(fileName, "a") do file
        println(file, "# PIMCID: $(Param.uid)\n")
        println(file, "# $invocation\n")
        
        println(file, "-------------- Begin Simulation Parameters ------------------------------------\n")
        println(file, "Command String\t\t\t:\t$invocation")
        println(file, "Ensemble\t\t\t:\tcanonical")
        println(file, "Simulation Type\t\t\t:\tPIMC")
        println(file, "Action Type\t\t\t:\tgsf")
        println(file, "Number of paths\t\t\t:\t$(Param.nPar)")
        println(file, "Interaction Potential\t\t:\tfree")   # Hard-coded for compatibility
        println(file, "External Potential\t\t:\tharmonic")  # Hard-coded for compatibility
        println(file, "Temperature\t\t\t:\t$(set["temp"])")
        println(file, "Chemical Potential\t\t:\t0.11000")   # Hard-coded for compatibility
        println(file, "Particle Mass\t\t\t:\t48.48")          # Hard-coded for compatibility
        println(file, "Number Time Slices\t\t:\t$(Param.nTsl)")
        println(file, "Specified Imaginary Time Step\t:\t$(set["tau"])")
        println(file, "Imaginary Time Step\t\t:\t$(set["tau"])")    # The actual code fudges here - it "rounds to the nearest nTsl integer and uses the corresponding tau from that
        println(file, "Imaginary Time Length\t\t:\t0.66667")# Hard-coded for compatibility
        println(file, "Initial Number Particles\t:\t$(Param.nPar)")   # Canonical, so fixed
        println(file, "Initial Density\t\t\t:\t0.10000")  # Hard-coded for (in)compatibility - is nPar/volume
        println(file, "Num. Broken World-lines\t\t:\t0")# Hard-coded for compatibility - not doing the World-Line stuff
        println(file, "Container Type\t\t\t:\tPrism") # Hard-coded for compatibility - definitely not touching this for a while
        println(file, "Container Dimensions\t\t:\t10.00000")    # Hard-coded for compatibility - this is where the dimension that the Production code was compiled with really counts
        println(file, "Container Volume\t\t:\t10.00000")    # See above - same story, but for volume
        println(file, "Lookup Table\t\t\t:\t1")   # I guess there's a lookup table used?
        println(file, "Maximum Winding Sector\t\t:\t1")
        println(file, "Initial Worm Constant\t\t:\t1.00000")    # Not using a worm
        println(file, "Initial CoM Delta\t\t:\t$(set["delta"])")
        println(file, "CoM Delta\t\t\t:\t$(set["delta"])")    # Because mine isn't changing for now
        println(file, "Bisection Parameter\t\t:\t16")   # Hard-coded; "m" in PIMC.jl's StagingMove()?
        println(file, "Update Length\t\t\t:\t256")    # Hard-coded for compatibility - I don't know what this affects :(
        println(file, "Potential Cutoff Length\t\t:\t10.00000") # Linear dimension size again
        println(file, "Bin Size\t\t\t:\t$(set["sweepsToBin"])")
        println(file, "Number EQ Steps\t\t\t:\t$(set["numEquil"])")
        println(file, "Number Bins Stored\t\t:\t$(set["numSamples"])")
        println(file, "Random Number Seed\t\t:\t12345") # Hard-coded - not actually the case
        println(file, "Virtual Window\t\t\t:\t5")

        println(file, "\n------------- End Simulation Parameters ---------------------------------------\n")


        println(file, "\n-------------- Begin Acceptance Data ------------------------------------------\n")

        println(file, "Total Rate\t\t\t:\t0.12345\n") # Bogus value

        println(file, "center of mass\t\t\t:\t$(Path.numAcceptCOM)\n")

        println(file, "bisection\t\t\t:\t$(Path.numAcceptStaging)\n")

        println(file, "\n-------------- End Acceptance Data --------------------------------------------\n")

        println(file, "\n-------------- Begin Estimator Data -------------------------------------------\n")
        # Figure out what's going on in this section
        println(file, "\n-------------- End Estimator Data ---------------------------------------------\n")
    end  
end

include("energy.jl")        # Contains the energy equations and estimators
include("action.jl")        # Contains the functions that compute action,etc.
include("PIMC.jl")          # The heart of the program, contains the methods
                            # that permute the particle(s)
include("expectations.jl")  # Contains the methods that compute <x>, <x^2>, etc.
#include("pyBinData.jl")     # Contains only the Python binning implementation

function GetCommandLineInvocation(ARGS)
    invocation = "julia $PROGRAM_FILE "
    for i = 1:length(ARGS)
        invocation *= ARGS[i] * " "
    end
    # There's technically that extra space at the end, but that's also in the production code's output, so...
    return invocation
end

# Optimize - https://docs.julialang.org/en/v1/manual/performance-tips
function main()
### Set-up for the simulation #################################################
    invocation = GetCommandLineInvocation(ARGS)

    set             = parse_commandline()
    tStamp          = Dates.value(Dates.now())
    lam             = 1/2 # A^2 * K -- hbar^2/(m k_B) -> m = hbar^2 / k_B
    temp            = set["temp"]       # Temperature to run simulation at
    numParticles    = set["numPart"]    # Number of particles
#    numTimeSlices   = set["numTime"]    # Number of beads along a line
    tau             = set["tau"]
    numEquilibSteps = set["numEquil"]   # Number of steps to equilibriate over
    observableSkip  = set["obsSkip"]    # Number of MC steps between observations
    numSamples      = set["numSamples"] # Defines numMCsteps 
    sweepsToBin     = set["sweepsToBin"]# Number of recorded MC sweeps to bin
    delta           = set["delta"]      # Sets the delta to be used for Shift()
    x_min           = set["Xmin"]       # Sets the position lower bound
    x_max           = set["Xmax"]       # Sets the position upper bound
    
    initDens        = numParticles / (x_max - x_min)    # Initial density

    # Number of MC steps to take in total
    numMCsteps = numEquilibSteps + observableSkip * numSamples
    
    # Imaginary time: beta / J (beta / M in Ceperley)
    numTimeSlices = trunc(Int64, 1/(tau*temp))
    #numMCtoBin  = 50
    # spatialBinWidth and numSpatialBins are for specific estimators; binSize is for binning the MC itself
    spatialBinWidth     = set["spatialBinWidth"]  # Width of bin for distribution binning in MC sampling
    numSpatialBins      = trunc(Int64, (x_max - x_min)/spatialBinWidth + 1)

    println("Simulation Parameters:")
    println("N\t\t= $(numParticles)")
    println("lambda\t\t= $(lam)")
    println("Temperature\t= $(temp)")
    println("numMCSteps\t= $(numMCsteps)\n")
    uid = uuid1()
    initDensity = 1.0
    imagTimeStep = 1.0
    
### cd to output directory ####################################################
    try
        cd("RESULTS/$tau")
    catch
        mkpath("RESULTS/$tau")
        cd("RESULTS/$tau")
    end

### Create files for output, write headers ####################################
#Production code naming scheme:$temp            $numParticles                       $initDensity            $imagTimeStep          $uid.dat"
#    file_name = "$(@sprintf("%06.3f",temp))-$(@sprintf("%04.0f",numParticles))-$(@sprintf("%06.3f",1.0))-$(@sprintf("%07.5f",1.0))-$uid.dat"
file_name = "$(@sprintf("%06.3f",temp))-$(@sprintf("%04.0f",numParticles))-$(@sprintf("%06.3f",initDens))-$(@sprintf("%07.5f",tau))-$uid.dat"
    
    println(file_name)

    estDatName = "ce-estimator-" * file_name
    WriteHeader(estDatName, "#PIMCID: $uid")
    WriteHeader(estDatName, "#\tE\t\t\tE2\t\t\tKE\t\t\tPE\t\t\tX\t\t\tX2")

### Initialize the main data structures and set random initial positions ######
    println("Initializing the data structures...")
    ke          = 0.0
    pe          = 0.0
    numAccCom   = 0
    numAccStag  = 0
    
    rng = MersenneTwister()
    
    beads = zeros(Float64, numTimeSlices, numParticles)
    # https://discourse.julialang.org/t/creating-an-array-matrix-with-a-specified-range-of-randomly-generated-numbers/15471
    for tSlice = 1:numTimeSlices
        for ptcl = 1:numParticles
            beads[tSlice,ptcl] = Shift(rng, delta)
        end
    end

    ### Setup the Paths object(s)
    determinants    = zeros(Float64, numTimeSlices, numParticles)
    potentials      = zeros(Float64, numTimeSlices, numParticles)
    numHistBins = 0

#TODO: I'VE NAMED THINGS STUPIDLY. UNSTUPIDIFY THEM. NAMELY, binWidth && numMCbins
    Prms = Params(numParticles,     #nPar
                    numTimeSlices,  #nTsl
                    lam,            #lam
                    tau,            #tau
                    x_min,          #left-most value
                    x_max,          #right-most value
                    delta,          #delta
                    numHistBins,    #numHistBins
                    sweepsToBin,    #sweepsToBin
                    spatialBinWidth,#binWidth
                    numSpatialBins, #numMCbins
                    numEquilibSteps,#numEquilibSteps
                    observableSkip, #observableSkip
                    numSamples,     #numSamples
                    file_name,      #baseName
                    uid)            #uuid
    Path = Paths(beads, determinants, potentials, ke, pe, numAccCom, numAccStag)

    # Write the log file so that there's at least that to work with
    logName = "ce-log-" * file_name
    WriteLogFileParameters(logName, Prms, Path, numMCsteps, set, invocation)

### Run the simulation proper #################################################
    println("Starting the simulation...")
    PIMC(Prms, Path, numMCsteps, set, rng)

### Collect and output any final results ######################################
    Path.numAcceptCOM       /= ( numEquilibSteps + sweepsToBin * numSamples )*numTimeSlices
    Path.numAcceptStaging /= ( numEquilibSteps + sweepsToBin * numSamples )*numTimeSlices
    logName = "ce-log-" * file_name
    WriteHeader(logName, "#PIMCID: $uid")
    WriteHeader(logName, "CoM Acceptance Ratio:\t\t$(Path.numAcceptCOM)")
    WriteHeader(logName, "Staging Acceptance Ratio:\t\t$(Path.numAcceptStaging)")
    println("Accepted CoM ratio: \t$(Path.numAcceptCOM)")
    println("Accepted Staging ratio: $(Path.numAcceptStaging)\n")
    
    cd("../..")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
