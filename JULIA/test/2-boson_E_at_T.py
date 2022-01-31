import numpy as np
import sys

def Single_Energy(j):
    return j + 0.5

def Partition(energy, T, L):
    sum = 0.0

    for i in range(L):
        sum += np.exp(-energy[i] / T)
    
    return sum

def main():
    T = float(sys.argv[1])
    
    hbar = 1.0
    w = 1.0
    kB = 1.0
    L = 10

    energy = np.zeros(L)
    for i in range(L):
        energy[i] = Single_Energy(i)

    Z1 = np.zeros(L)
    Z2 = np.zeros(L)
    eZ1 = np.zeros(L)
    eZ2 = np.zeros(L)
    for i in range(L):
        for j in range(L):
            #Z1[i] = Partition(energy, T, L)
            Z1[i] += np.exp(-(energy[i] + energy[j])/T)
            eZ1[i] += (energy[i] + energy[j]) * np.exp(-(energy[i] + energy[j])/T)
        #Z2[i] = Partition(2.0*energy, T, L)
        Z2[i] = np.exp(-2.0 * energy[i] / T)
        eZ2[i] = 2.0 * energy[i] * np.exp(-2.0 * energy[i] / T)
    
    denom = sum(Z1) + sum(Z2)
    
    numer = 0.0
    for i in range(L):
        for j in range(L):
            numer += (energy[i] + energy[j]) * np.exp(-(energy[i] + energy[j])/T)

        numer += 2.0 * energy[i] * np.exp(-2 * energy[i] / T)


    E = numer / denom

    print(f"E is approximately: {E}")

if __name__ == "__main__":
    main()
