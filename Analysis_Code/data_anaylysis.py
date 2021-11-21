# For a given set of paramters, combines the corresponding data files into one 
# file

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mpm
import sys as sys

def ComputeSampleAverage(x):
    average = np.zeros(len(x[0]))
    N = len(x)
    
    for i in range(len(x[0])):
        for j in range(N):
            average[i] += x[j][i]

    return average/N

# Computes the Jackknife average of each element
def ComputeJackAverage(x):
    jAverage = np.zeros_like(x)
    N = len(x)
    print(N)
    for j in range(len(x[0])):
        for i in range(N):
            for k in range(N):
                jAverage[i][j] += x[k][j]
                if k == i:
                    jAverage[i][j] -= x[k][j]

    return jAverage/(N-1)

# Compute the Jackknife variance

def ComputePseudoValue(average, jAverage):
    pseu = np.zeros_like(jAverage)
    N = len(jAverage)

    for i in range(len(jAverage)):
        for j in range(len(jAverage[0])):
            pseu[i][j] = average[j] + (N - 1)*(average[j] - jAverage[i][j])

    return pseu

def ComputePseudoAverage(pseu):
    psAverage = np.zeros_like(pseu[0])
    
    for i in range(len(pseu[0])):
        for j in range(len(pseu)):
            psAverage[i] += pseu[j][i]

    return psAverage/(len(pseu))

def ComputePseudoVariance(pseu,psAverage):
    psVariance = np.zeros_like(psAverage)
    N = len(pseu)

    for i in range(len(pseu[0])):
        for j in range(N):
            psVariance[i] += pow( pseu[j][i] - psAverage[i] ,2)
    # Something about dividing by the square root of N?
    return psVariance/np.sqrt((N - 1))

def main():
    temp = 1.50
    steps = 3000
    numBeads = 50
    numPart = 1
    obs = 50

    ObservationSteps = [50,75,100,125,150]
    Temps = [.5,.75,1.0,1.25,1.5]
    NumberBeads = [25,50,75,100,125]
    MCSteps = [3000,4000,5000,6000,7000]

    x = np.zeros((100,5))
    
    # Iterates over the files by the desired aspect
    for step in MCSteps:
        name = "data_T_" + str(temp) + "-Eq_" + str(step) + "-Obs_" + str(obs) + "-nB_" + str(numBeads) + "-nP_1.txt"
        print(f"Opening file {name} ...")
        with open(name,"r") as file1:
            # Read the data into the dat array - it will be used to make plots and do the analysis for later
            dat = file1.read().splitlines()
        file1.close()

        for j in range(3,103):
            x[j-3][MCSteps.index(step)] = float(dat[j])
        
    exactE = 1/2*mpm.tanh(1/(2*temp))
    print(f"The exact energy is {exactE}")
    # phi_n(X)
    average = ComputeSampleAverage(x)
    print(f"The average is {average}")
    # phi_(n-1)(X)_i
    jAverage = ComputeJackAverage(x)
    #print(jAverage)
    pseu = ComputePseudoValue(average,jAverage)
    #print(pseu)
    psAverage = ComputePseudoAverage(pseu)
    print(f"The \"pseudo-average\" is {psAverage}")
    psVariance = ComputePseudoVariance(pseu,psAverage)
    print(f"The \"pseudo-variance\" is {psVariance}")

################################################################################

    # Plot the stuff
    fig = plt.figure()
    plt.scatter(MCSteps, psAverage)
    plt.errorbar(MCSteps, psAverage, yerr=psVariance,ls='none')
    plt.xlabel("# of MC Steps")
    plt.ylabel("Average Energy ()")
    
    plt.hlines(float(exactE),MCSteps[0],MCSteps[len(MCSteps)-1],linestyle="dashed",label=r'$E_{exact}$')

    figName = "Energy_" +str(temp)+ "-nB_" + str(numBeads) + "-obs_" + str(obs) + ".png"
    fig.savefig(figName)

    plt.show()

if __name__ == "__main__":
    main()
