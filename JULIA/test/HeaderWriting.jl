# Writes the header of the given file, filname, with the given string
# Technically, it only writes a single string at a time, so it can be used for
# writing single lines.
function WriteHeiader(fileName::String, string::String)
    open(fileName, "a") do file
        println(file, string)
    end
end

# Writes out a whole matrix of data to a given file
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
