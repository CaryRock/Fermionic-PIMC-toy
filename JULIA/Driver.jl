# Notes
# * The Max Graves thesis has the derivation of the virial implementation of the
# KE estimator - unnecessary for this code, but will be interesting/useful later

using Dates
using Random
using ArgParse
using Printf
using ProgressBars  # Entirely vain - just gives the pretty progress bar
using UUIDs         # I sort-of despise this language. I really kind of do.
struct Params
    nPar::Int64             # Number of particles
    nTsl::Int64             # Number of time slices
    lam::Float64            # System parameters - hbar^2/(2mk_B) = 1/2
    tau::Float64            # beta/J (beta/M in Ceperley-ese)
    x_a::Float64            # Arbitrary left edge
    x_b::Float64            # Arbitrary right edge
    delta::Float64          # Width of possible shift for a bead
    numHistBins::Int64      # Number of bins for histogram
    numBins::Int64
    # The following two should be specific to relevant estimators, not global
    binWidth::Float64       # Width of bins for binning the distribution
    numMCbins::Int64        # Number of bins for MC binning
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
    numAcceptCOM::Int64
    numAcceptStaging::Int64
end

function parse_commandline()
    set = ArgParseSettings()

    @add_arg_table set begin
### Required values ###########################################################
        "--temp", "-T"
            help = "Temperature in K"
            required = true
            arg_type = Float64
        "--numPart", "-N"
            help = "Number of particles"
            required = true
            arg_type = Int64
        "--numTime", "-J"
            help = "Number of time slices"
            required = true
            arg_type = Int64
        "--numEquil", "-E"
            help = "Number of steps to equilibriate over"
            required = true
            arg_type = Int64
        "--obsSkip", "-O"
            help = "Number of steps between observations"
            required = true
            arg_type = Int64
        "--numSamples", "-S"
            help = "Number of desired samples to take"
            required = true
            arg_type = Int64
            default = 201
### Options for Setting Variables #############################################
        "--binWidth", "-b"
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

### Optional Output Flags #####################################################
# These currently don't do anything, but may be reimplemented later
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
@inline function CutOff(beadA::Float64, beadB::Float64, cutoff::Float64=10^-12)
    return ( abs(beadA - beadB) < cutoff )
end

# Controls the width of the possible shift in position for a bead
@inline function Shift(delta::Float64  = 1.0)
    return delta * (-1.0 + 2.0 * rand(MersenneTwister()))
end

# Writes the header of the given file, filename, with the given string
# Technically, it only writes a single string at a time, so it can be used for 
# writing single lines
function WriteHeader(filename::String, string)
    open(filename, "a") do file
        println(file, string)
    end
end

# Writes out a whole matrix work of data (data) to a given file (filename)
function WriteOutMat(filename::String, data, xLim::Int64 = length(data), yLim::Int64 = length(data[1]))
    open(filename, "a") do file
        for t = 1:xLim
            for p = 1:yLim
                println(file, data[t,p])
            end
        end
    end
end

include("energy.jl")        # Contains the energy equations and estimators
include("action.jl")        # Contains the functions that compute action,etc.
include("PIMC.jl")          # The heart of the program, contains the methods
                            # that permute the particle(s)
include("expectations.jl")  # Contains the methods that compute <x>, <x^2>, etc.
#include("pyBinData.jl")     # Contains only the Python binning implementation

# Optimize - https://docs.julialang.org/en/v1/manual/performance-tips
function main()
### Set-up for the simulation #################################################
    set             = parse_commandline()
    tStamp          = Dates.value(Dates.now())
    lam             = 1/2 # A^2 * K -- hbar^2/(m k_B) -> m = hbar^2 / k_B
    temp            = set["temp"]       # Temperature to run simulation at
    numParticles    = set["numPart"]    # Number of particles
    numTimeSlices   = set["numTime"]    # Number of beads along a line
    numEquilibSteps = set["numEquil"]   # Number of steps to equilibriate over
    observableSkip  = set["obsSkip"]    # Number of MC steps between observations
    numSamples      = set["numSamples"] # Defines numMCsteps 
    delta           = set["delta"]      # Sets the delta to be used for Shift()
    x_min           = set["Xmin"]       # Sets the position lower bound
    x_max           = set["Xmax"]       # Sets the position upper bound

    # Number of MC steps to take in total
    # TODO: IMPLEMENT A while() LOOP SO THAT THIS IS UNNECESSARY
    numMCsteps = numEquilibSteps + observableSkip * numSamples
    
    # Imaginary time: beta / J (beta / M in Ceperley)
    tau = 1/(numTimeSlices * temp)

    #numMCtoBin  = 50
    # binWidth and numMCbins are for specific estimators; binSize is for binning the MC itself
    binWidth    = set["binWidth"]  # Width of bin for distribution binning in MC sampling
    numMCbins   = trunc(Int64, (x_max - x_min)/binWidth)

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
        cd("RESULTS/$(set["numTime"])")
    catch
        mkdir("RESULTS/$(set["numTime"])")
        cd("RESULTS/$(set["numTime"])")
    end

### Create files for output, write headers ####################################
#Production code naming scheme:$temp            $numParticles                       $initDensity            $imagTimeStep          $uid.dat"
#    file_name = "$(@sprintf("%06.3f",temp))-$(@sprintf("%04.0f",numParticles))-$(@sprintf("%06.3f",1.0))-$(@sprintf("%07.5f",1.0))-$uid.dat"
    file_name = "$(@sprintf("%06.3f",temp))-$(@sprintf("%04.0f",numParticles))-$(@sprintf("%04.0f",numTimeSlices))-$uid.dat"
    
    println(file_name)
#    file_name = "data_T_$temp-Eq_$numEquilibSteps-Obs_$observableSkip-nB_$numTimeSlices-nP_$numParticles-$tStamp.dat"

    estDatName = "ce-estimator-" * file_name
    WriteHeader(estDatName, "#PIMCID: $uid")
    WriteHeader(estDatName, "#\tE\t\t\tE2\t\t\tKE\t\t\tPE\t\t\tX\t\t\tX2")

### Initialize the main data structures and set random initial positions ######
    println("Initializing the data structures...")
    ke          = 0.0
    pe          = 0.0
    numAccCom   = 0
    numAccStag  = 0
    beads = zeros(Float64, numTimeSlices, numParticles)
    # https://discourse.julialang.org/t/creating-an-array-matrix-with-a-specified-range-of-randomly-generated-numbers/15471
    for tSlice = 1:numTimeSlices
        for ptcl = 1:numParticles
            beads[tSlice,ptcl] = Shift(delta)
        end
    end

    ### Setup the Paths object(s)
    determinants    = zeros(Float64, numTimeSlices, numParticles)
    potentials      = zeros(Float64, numTimeSlices, numParticles)
    numHistBins = 0
    numBins = 50

#TODO: I'VE NAMED THINGS STUPIDLY. UNSTUPIDIFY THEM. NAMELY, binWidth && numMCbins
    Prms = Params(numParticles,     #nPar
                    numTimeSlices,  #nTsl
                    lam,            #lam
                    tau,            #tau
                    x_min,          #x_a
                    x_max,          #x_b
                    delta,          #delta
                    numHistBins,    #numHistBins
                    numBins,        #numBins
                    binWidth,       #binWidth
                    numMCbins,      #numMCbins
                    numEquilibSteps,#numEquilibSteps
                    observableSkip, #observableSkip
                    numSamples,     #numSamples
                    file_name,      #baseName
                    uid)            #uuid
    Path = Paths(beads, determinants, potentials, ke, pe, numAccCom, numAccStag)

### Run the simulation proper #################################################
    println("Starting the simulation...")
    PIMC(Prms, Path, numMCsteps, set)

### Collect and output any final results ######################################
    println("Accepted CoM moves: $(Path.numAcceptCOM/numMCsteps)")
    println("Accepted Staging moves: $(Path.numAcceptStaging/numMCsteps)\n")
    
    cd("../..")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
