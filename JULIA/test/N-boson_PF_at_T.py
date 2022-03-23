import numpy as np
import sys

class PF:
    # Partition function class. Each should have at least a Z and an eZ
    # for recursion
    def __init__(self, L, Potential, Interaction):
        self.Z = np.zeros(L)
        self.eZ = np.zeros(L)
        self.PotentialE = Potential
        self.InteractionE = Interaction

def Single_Energy(j):
    return j + 0.5

def Partition(energy, T, L):
    sum = 0.0

    for i in range(L):
        sum += np.exp(-energy[i] / T)
    
    return sum

def main():
    T = float(sys.argv[1])
    N = int(sys.argv[2])

    hbar = 1.0
    w = 1.0
    kB = 1.0
    L = int(20.0 * T)

    energy = np.zeros(L)
    for i in range(L):
        energy[i] = Single_Energy(i)
    
    #Zarray = [PF(L) for i in range(N)]
    
    Z1 = np.zeros(L)    # Zarray[0].Z
    Z2 = np.zeros(L)    # Zarray[1].Z
    eZ1 = np.zeros(L)   # Zarray[0].eZ
    eZ2 = np.zeros(L)   # Zarray[1].eZ
    for i in range(L):
        for j in range(L):
            Z1[i] += np.exp(-(energy[i] + energy[j])/T)
            eZ1[i] += (energy[i] + energy[j]) * np.exp(-(energy[i] + energy[j])/T)
        Z2[i] = np.exp(-2.0 * energy[i] / T)
        eZ2[i] = 2.0 * energy[i] * np.exp(-2.0 * energy[i] / T)
    
    denom = sum(Z1) + sum(Z2)
    numer = 0.0
    for i in range(L):
        for j in range(L):
            numer += (energy[i] + energy[j]) * np.exp(-(energy[i] + energy[j])/T)

        numer += 2.0 * energy[i] * np.exp(-2 * energy[i] / T)


    E = numer / denom

    print(f"{T}\t{E}")

if __name__ == "__main__":
    main()
