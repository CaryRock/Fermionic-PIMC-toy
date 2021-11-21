using PyCall

function PyBinData(energyTrace::Array{Float64},kineticTrace::Array{Float64}, 
            potentialTrace::Array{Float64}, binSize::Int64, numParticles::Int64,
            temp::Float64,equil::Int64,obs::Int64, numTimeSlices::Int64, 
            numProcesses::Int64, file_name::String)
    py"""
    import numpy as np
    #import mpmath as mp
    pyTemp = $temp
    pyEquil = $equil
    pyObs = $obs
    pyNumTimeSlices = $numTimeSlices
    pyNumParticles = $numParticles
    Energy = $vec($energyTrace)
    Kinetic = $vec($kineticTrace)
    Potential = $vec($potentialTrace)
    pyBinSize = $binSize
    numBins = int(1.0 * len(Energy)/pyBinSize)
    slices = np.linspace(0,len(Energy),numBins + 1,dtype=int)
    binnedEnergy = np.add.reduceat(Energy,slices[:-1]) / np.diff(slices)
    sumEnergy = np.mean(binnedEnergy) #/11605
    errEnergy = np.std(binnedEnergy)/np.sqrt(numBins-1) #/11605
    # 11605 K per 1 eV
    #print('Energy = %8.4f +/- %6.4f' % (sumEnergy, errEnergy))
    #print('Eexact = 1/2*coth(1/( 2T )) = ', $numParticles*1/4*mp.coth(1/(2 * $temp)))
    name = $file_name
    with open(name,"a") as output:
        print(sumEnergy,file=output)
    output.close()
    """
end
