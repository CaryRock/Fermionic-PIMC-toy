#! /usr/bin/env python

def main():
    # sample while() loop
    sampleCount     = 0     # Counts how many samples taken
    binCount        = 0     # Counts how many times gone through binning
    numSamples      = 5     # The desired number of samples
    equilibSkip     = 0
    steps           = equilibSkip + 1 # The current MC "sweep" being performed
    observableSkip  = 1     # This requires -O 1 at the least
    nTsl            = 3     # Number of beads on a line
    nPar            = 1     # Number of particles (lines)
    binSize         = 5     # Number of samples to bin together before writing out
    finalX1         = 0     # Write out value for X1
    finalX2         = 0     # Write out value for X2
    x1_ave = 0
    x2_ave = 0

    while sampleCount < numSamples:
        # Update(Param, Path)
    
        # If observableSkip steps have occured, record information
        if (steps % observableSkip == 0):
            binCount += 1
    
            for i in range(nTsl):
                for j in range(nPar):
                    x1_ave += -1 + i
                    x2_ave = (-1 + i) * (-1 + i)

            if (binCount % binSize == 0):
                sampleCount += 1
                x1_ave /= binSize
                x2_ave /= binSize
                # Write out x1_ave, x2_ave, etc.
                print(f"Sample #: {sampleCount}\tbinCount #: {binCount}")
                print(f"X1_ave: {x1_ave}\tX2_ave: {x2_ave}")
                print("")

                # After recording, zero out the accumulators as well as binCount (though it doesn't technically need to be)
                x1_ave      = 0
                x2_ave      = 0
                binCount    = 0

if __name__ == "__main__":
    main()
