import numpy as np
import sys

def Single_Energy(j):
    return j + 0.5

# Partition() through Deriv_SofK() can be condensed
def Partition(energy, T, L):
    sum = 0.0

    for i in range(L):
        sum += np.exp(-energy[i] / T)
    
    return sum

def SofK(energy, T, L, mult):
    sum = 0.0

    for i in range(L):
        sum += np.exp(-energy[i] * mult / T)

    return sum

def Deriv_Partition(energy, T, L):
    sum = 0.0

    for i in range(L):
        sum += energy[i] * np.exp(-energy[i] / T)

    return -sum

def Deriv_SofK(energy, T, L, mult):
   sum = 0.0

   for i in rage(L):
       sum += energy[i] * np.exp(-energy[i] * mult / T)

    return -mult * sum

def main():
    T = float(sys.argv[1])
    N = int(sys.argv[2])

    hbar = 1.0
    w = 1.0
    kB = 1.0
    L = 10  # Cutoff - the number of energy levels to compute to approximate "infinity"

    energy = np.zeros(L)
    for i in range(L):
        energy[i] = Single_Energy(i)
    
    Zarray = [0.0 for i in range(N + 1)] # Array to store PFs - N+1 because 0 -> N









    
    Z1 = np.zeros(L)    # Zarray[0].Z
    Z2 = np.zeros(L)    # Zarray[1].Z
    eZ1 = np.zeros(L)   # Zarray[0].eZ
    eZ2 = np.zeros(L)   # Zarray[1].eZ
    for i in range(L):
        for j in range(L):
            #Zarray[i].Z += np.exp(-(energy[i] + energy[j])/T)
            Z1[i] += np.exp(-(energy[i] + energy[j])/T)
            #Zarray[i].eZ+= (energy[i] + energy[j]) * np.exp(-(energy[i] + energy[j])/T)
            eZ1[i] += (energy[i] + energy[j]) * np.exp(-(energy[i] + energy[j])/T)
        #Z2[i] = Partition(2.0*energy, T, L)
        Z2[i] = np.exp(-2.0 * energy[i] / T)
        eZ2[i] = 2.0 * energy[i] * np.exp(-2.0 * energy[i] / T)
    
    denom = sum(Z1) - sum(Z2)
    # denon = np.sum(Zarray)
    numer = 0.0
    for i in range(L):
        for j in range(L):
            numer += (energy[i] + energy[j]) * np.exp(-(energy[i] + energy[j])/T)

        numer -= 2.0 * energy[i] * np.exp(-2 * energy[i] / T)


    E = numer / denom

    print(f"E is approximately: {E}")

if __name__ == "__main__":
    main()
