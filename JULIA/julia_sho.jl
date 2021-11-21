# Notes
# * Implement the "Levy Bridge" (Staging Move) code to help with convergence - 
# it can basically be lifted from the sho-pimc.py code
# * The Ceperly paper has the derivation for the N*k_B*T*J term appearing in the 
# KE estimator
# * The Max Graves thesis has the derivation of the Virial implementation of the
# KE Estimator - unnecessary for this code, but will be interesting

using Dates
using Random
using ArgParse

ver="11"

mutable struct Paths
    numParticles::Int64
    numTimeSlices::Int64
    beads::Array{Float64,}
    lam::Float64
    tau::Array{Float64}
    determinants::Array{Float64}
    potentials::Array{Float64}
    KE::Float64
    PE::Float64
    numAcceptCOM::Int64
    numAcceptStaging::Int64
end

#struct Params
#    nPar::Int64        # Number of particles
#    nTsl::Int64        # Number of time slices
#    lam::Float64       # System parameters - hbar^2/(2mk_B) = 1/2
#    tau::Array{Float64}# beta/J
#    x_a::Float64       # Arbitrary left edge
#    x_b::Float64       # Arbitrary right edge
#    numBins::Int64     # Number of buckets for histogram?
#    numMCbins::Int64   # Number of bins for MC binning
#end

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
            arg_type = Int
        "--numTime", "-J"
            help = "Number of time slices"
            required = true
            arg_type = Int
        "--numEquil", "-E"
            help = "Number of steps to equilibriate over"
            required = true
            arg_type = Int
        "--obsSkip", "-O"
            help = "Number of steps between observations"
            required = true
            arg_type = Int

### Options for Setting Variables #############################################
        "--binWidth", "-b"
            help = "Sets the desired width of bins"
            arg_type = Float64
            default = 0.01
        "--delta", "-d"
            help = "Sets the delta for the CoM move"
            arg_type = Float64
            default = 1.0
        "--Xmin"
            help = "Sets the minimum value for X - can be negative"
            arg_type = Float64
            default = -10.0
        "--Xmax"
            help = "Sets the maximum value for X - can be negative"
            arg_type = Float64
            default = 10.0
### Optional Output Flags #####################################################
#=
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
=#
    end

    return parse_args(set)
end

# The external potential to use
@inline function ExtPotential(lam::Float64,R::Float64)
    return 1/(4 * lam ) * R * R
end

# Computes the modular time slice to use
# +1 because Julia is 1-indexed, not 0-indexed like Python or C
@inline function ModTslice(tSlice::Int64,numSlices::Int64)
    return (( abs(tSlice) - 1 + numSlices ) % numSlices ) + 1
end

# Flags whether or not to truncate the calculation
@inline function CutOff(beadA::Float64,beadB::Float64,cutoff=10^-12)
    return (abs(beadA - beadB) < cutoff)
end

@inline function Shift(delta=0.75)
    return delta*(-1.0 + 2.0 * rand(MersenneTwister()))
end

function WriteHeader(filename, string)
    open(filename, "a") do file
        println(file, string)
    end
end

function WriteOutMat(filename, data, xLim=length(data), yLim=length(data[1]))
    open(filename, "a") do file
        for t = 1:xLim
            for p = 1:yLim
                println(file, data[t,p])
            end
        end
    end
end

include("energy.jl")        # Contains the energy equations
include("action.jl")        # Contains the functions that compute action, etc.
include("PIMC.jl")          # The heart of the program; contains the methods 
                            # that permute the particle(s)
include("expectations.jl")  # Contains the methods that compute <x>, <x^2>, etc.
include("pyBinData.jl")     # Contains simply Python binning

# Optimize - https://docs.julialang.org/en/v1/manual/performance-tips/
function main()
    set             = parse_commandline()
    tStamp          = Dates.value(Dates.now())
    lam             = 1/2 # A^2 * K - \hbar^2/m k_b -> m = hbar^2 / k_B
    initialT        = 0.01                      # Starting temperature, K
    finalT          = set["temp"]       #parse.(Float64, ARGS[1])  # Final temperature, K
    numParticles    = set["numPart"]    #parse.(Int64, ARGS[2])    # Number of particles    
    numTimeSlices   = set["numTime"]    #parse.(Int64, ARGS[3])    # Number of beads along a line
    numEquilibSteps = set["numEquil"]   #parse.(Int64, ARGS[4])
    observableSkip  = set["obsSkip"]    #parse.(Int64, ARGS[5])
    multiple        = 11                # Defines numMCSteps
    delta           = set["delta"]      # Sets delta for Shift() function
    numMCSteps = multiple*numEquilibSteps + observableSkip 
                                                # Number of MC steps to take
    leng = 1 + trunc(Int, (multiple - 1) * numEquilibSteps / observableSkip)
        # How long to make the various recording data arrays

    temp = finalT
    # Imaginary time: \beta/M
    tau = 1/(numTimeSlices * finalT)

    binSize = 50        # Data binning in Python
    binWidth= 0.01      # MC binning

    println("Simulation Parameters:")
    println("N\t\t= $(numParticles)")
    println("lambda\t\t= $(lam)")
    println("initialT\t= $(initialT)")
    println("finalT\t\t= $(finalT)")
    println("numMCSteps\t= $(numMCSteps)")
    println()

    file_name = "data_T_$finalT-Eq_$numEquilibSteps-Obs_$observableSkip-nB_$numTimeSlices-nP_$numParticles-$tStamp.dat"
    #WriteOutMat(file_name, [numEquilibSteps, observableSkip, finalT])

    enDatName = "raw_energy_" * file_name
    WriteHeader(enDatName, "#\tE\t\t\tKE\t\t\tPE")
    
    energyTrace     = zeros(Float64, leng)
    potentialTrace  = zeros(Float64, leng)
    kineticTrace    = zeros(Float64, leng)
    x1_ave          = 0.0
    x2_ave          = 0.0

    println("Starting run $iter of $numProcesses")

### Initialize the main data structure and set random initial positions #######
    beads = zeros(Float64,numTimeSlices,numParticles)
    # https://discourse.julialang.org/t/creating-an-array-matrix-with-a-specified-range-of-randomly-generated-numbers/15471
    for tSlice = 1:numTimeSlices
        for ptcl = 1:numParticles
             beads[tSlice,ptcl] = Shift(delta)   #0.5 * (-1.0 + 2.0 * rand())
        end
    end

    ### Setup the Paths
    determinants = zeros(Float64,numTimeSlices,numParticles)
    potentials = zeros(Float64,numTimeSlices,numParticles)

    Path = Paths(numParticles, numTimeSlices, beads, lam, tau, determinants,
        potentials, 0.0, 0.0, 0, 0)
    
# TODO: Basically rewrite this whole file
    energyTrace, kineticTrace, potentialTrace = PIMC(numMCSteps, Path, iter, numEquilibSteps, observableSkip, file_name, multiple, binWidth, set)
    
### Write out raw energy data #################################################
    # Write raw energy data out to file
    enDatName = "raw_energy_" * file_name
    open(enDatName, "a") do file
        for dat=1:length(energyTrace)
            println(file, energyTrace[dat], "\t", kineticTrace[dat], "\t",
                    potentialTrace[dat])
        end
    end
    
### data and end program ##################################################

    # Copied Python binning from SHO
#=
    PyBinData(energyTrace, kineticTrace, potentialTrace, binSize,
              numParticles, temp[iter], numEquilibSteps, observableSkip,
              numTimeSlices, file_name)
=#

    println("Accepted CoM: $(Path.numAcceptCOM/numMCSteps)")
    println("Accepted Staging: $(Path.numAcceptStaging/numMCSteps)\n")
       
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
