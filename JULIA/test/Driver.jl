# Notes
# * The Max Graves thesis has the derivation of the virial implementation of the
# KE estimator - unnecessary for this code, but will be interesting/useful later
using MKL           # Intel's linear algebra library - fastfast
using LinearAlgebra # Linear algebra library - functions, etc.
using Dates         # To get the date
using Random        # Library of pRNGs and related
using ArgParse      # To parse commandline arguments
using Printf        # To help write files
using ProgressBars  # Entirely vain - just gives the pretty progress bar
using UUIDs         # I sort-of despise this language. I really kind of do.

#Using Threads       # Testing - not used
#using TimerOutputs  # Optimization profiling

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
    detMat::Array{Float64, 3}   # N x N particles, nTsl slices; Julia's built different
    determinants::Array{Float64}# Holds evaluated abs(ln(array))
    potentials::Array{Float64}  # Computed Potentials contribution to action
    KE::Float64                 # Kinetic Energy
    PE::Float64                 # Potential Energy
    numAcceptCOM::Float64       # Number of Accepted from COM operation
    numAcceptStaging::Float64   # Number of Accepted from Staging operation
end

#=
Parses commandline arguments, returns set of argumnets to main program for 
simulation use.
=#
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

#= 
The external potential to use
Input: coefficient, location
If more start being added, this should go into a separate file with the other
potentials
=#
@inline function ExtPotential(lam::Float64, R::Float64)
    return 1/(4 * lam) * R * R
end

#=
# Computes the modular time slice to use
# +1 because Julia is 1-indexed, not 0-indexed like Python or C
Input: value to mod, base of ring
=#
@inline function ModTslice(tSlice::Int64, numSlices::Int64)
    return (( abs(tSlice) - 1 + numSlices ) % numSlices ) + 1
end

#=
# Flags whether or not to truncate the calculation
Input: Two locations, minimum distance before cutoff
=#
@inline function CutOff(beadA::Float64, beadB::Float64, cutoff::Float64=10^-4)
    return ( abs(beadA - beadB) < cutoff )
end

#=
# Controls the width of and shift in the position for a bead
Input: Desired pRNG, scaling factor for shift
=#
@inline function Shift(rng::MersenneTwister, delta::Float64  = 1.0)
    return delta * (-1.0 + 2.0 * rand(rng))
end

### Program files to include ##################################################
include("HeaderWriting.jl") # Contains all the file writing functions
include("energy.jl")        # Contains the energy equations and estimators
include("action.jl")        # Contains the functions that compute action,etc.
include("PIMC.jl")          # The heart of the program, contains the methods
                            # that permute the particle(s)
include("expectations.jl")  # Contains the methods that compute <x>, <x^2>, etc.
#include("pyBinData.jl")     # Contains only the Python binning implementation

#=
Gets the commandline invocation for later writing to a file.
Input: Commandline
=#
function GetCommandLineInvocation(ARGS)
    invocation = "julia $PROGRAM_FILE "
    for i = 1:length(ARGS)
        invocation *= ARGS[i] * " "
    end
    # There's technically that extra space at the end, but that's also in the
    # production code's output, so...
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
    #numTimeSlices   = set["numTime"]    # Number of beads along a line
    tau             = set["tau"]        # Distance between imaginary time steps
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
    initDensity = 1.0   # For file writing to match production code's - unused
    imagTimeStep = 1.0  # For file writing to match production code's - unused
    
### cd to output directory ####################################################
    try
        cd("RESULTS/$tau")
    catch
        mkpath("RESULTS/$tau")
        cd("RESULTS/$tau")
    end

### Create files for output, write headers ####################################
    #Production code naming scheme:$temp            $numParticles                       $initDensity            $imagTimeStep          $uid.dat"
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
    # Tensors go Row x Column x "slice"
    detMat            = zeros(Float64, numParticles, numParticles, numTimeSlices)
    #@printf "Size of detMat tensor: %s\n" string(size(detMat))
    determinants    = zeros(Float64, numTimeSlices)
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
    Path = Paths(beads, detMat, determinants, potentials, ke, pe, numAccCom, numAccStag)

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

### Program run ###############################################################
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end