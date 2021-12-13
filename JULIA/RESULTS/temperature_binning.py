#! /usr/bin/env python3

import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import mpmath as mpm

def Error_Binning(data):
    nBins = data.size
    delta = np.std(data)/np.sqrt(nBins)

    start = nBins % 2
    binned_data = 0.5 * ( data[start::2] + data[start + 1::2] )
    
    return delta, binned_data

def Compute_Exact(tempRange):
    return 1/2 * mpm.coth(1/(2*tempRange))

def main():
    L = int(sys.argv[1])
#    e = int(sys.argv[2])

    temps = np.array([0.3,0.5,0.75,1.0,1.25,1.5])
    fig1 = plt.figure(1)
    fig2 = plt.figure(2)

    files = sorted(glob.glob("ce-estimator-*"))
    nFiles = len(files) 
    
    tempRange = np.arange(0.15,1.65, 0.025)
    exact = np.zeros_like(tempRange)

    for i in range(len(tempRange)):
        exact[i] = Compute_Exact(tempRange[i])
    
    plt.figure(1)
    plt.plot(tempRange, exact, 'r--', label="Exact")
    
    plt.figure(2)
    plt.plot(tempRange, exact, 'r--', label="Exact")

    nSamples = 0
    for i in range(nFiles):
        with open(files[i], "r") as f:
            lines = (line for line in f if not line.startswith('#'))
            data = np.loadtxt(lines)
            # 'data' is simply the array of data in the file - E E2 KE PE X X2
            # Since this is over temperature, will be binning E over temperature
            
            energyData  = data[:,0]
            x2Data      = data[:,5]
            nSamples    = len(energyData)

            nLevels     = int(np.log2(energyData.size))

            for j in range(nLevels):
                errorBars, energyData   = Error_Binning(energyData)
                x2Bars, x2Data          = Error_Binning(x2Data)

            plt.figure(1)
            plt.errorbar(temps[i], energyData, yerr=errorBars, fmt='.', label=f"T_data = {temps[i]}")
            
            plt.figure(2)
            plt.errorbar(temps[i], x2Data, yerr=x2Bars, fmt='.', label=f"T_data = {temps[i]}")

        f.close()

    plt.figure(1)
    plt.title(f"Energy vs. Temperature @ M = {L}; # Samples = {nSamples}")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Energy (K)")
    plt.legend()

    plt.figure(2)
    plt.title(f"X2 vs. Temperature @ M = {L}; # Samples = {nSamples}")
    plt.xlabel("Temperature (K)")
    plt.ylabel("X2 (A^2)")
    plt.legend()

    plt.show()

if __name__ == "__main__":
    main()
