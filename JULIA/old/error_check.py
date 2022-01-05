import numpy as np
import matplotlib.pyplot as plt
import sys
import mpmath as mpm

def Exact(temp):
    return 1/2*mpm.coth(1/(2*temp))

def main():
    energyFile = sys.argv[1]
    errFile = sys.argv[2]
    
    dataEnergy = []
    dataErr = []
    with open(energyFile,"r") as f:
        for line in f.readlines():
            dataEnergy.append(float(line))
    f.close()

    with open(errFile,"r") as g:
        for line in g.readlines():
            dataErr.append(float(line))
    g.close()

    _energies = dataEnergy
    _errs = dataErr
    
    x = np.array([0.1, 0.5, 1.0, 5.0])
    N = len(x)
    y = np.zeros(N)
    errs = np.zeros(N)

    for i in range(N):
        y[i] = float(_energies[i])
        errs[i] = float(_errs[i])

    t = np.linspace(0.1,5.0,50)
    z = np.zeros(len(t))

    for j in range(len(t)):
        z[j] = Exact(t[j])

    fig = plt.figure()
    plt.plot(x, y, 'bo')
    plt.errorbar(x, y, yerr=errs,ls='none')
    plt.plot(t, z, 'r--')
    plt.xlabel("System Temperature (K)")
    plt.ylabel("Energies (K)")
    
    plt.legend(["Computed Data","Exact Computation"])
    fig.savefig("present.png")
    
    plt.show()

if __name__ == "__main__":
    main()
