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

    for i in range(L):
        sum += energy[i] * np.exp(-energy[i] * mult / T)

    return (-mult * sum)

def main():
    T = float(sys.argv[1])
    N = int(sys.argv[2])
    Nplus = N + 1

    hbar = 1.0
    w = 1.0
    kB = 1.0
    L = 10  # Cutoff - the number of energy levels to compute to approximate "infinity"
    eta = -1
    energy = np.zeros(L)
    for i in range(L):
        energy[i] = Single_Energy(i)
    
    Z_array = [1.0 for i in range(Nplus)] # Array to store PFs - N+1 because 0 -> N
    dZdB_array = [1.0 for i in range(Nplus)]
    Sk_array = [1.0 for i in range(Nplus)]
    dSkdB_array = [1.0 for i in range(Nplus)]

    for i in range(Nplus):
        Z_array[i] = Partition((N - i) * energy, T, L)   # Produces Zarray[Z_N, Z_N-1, ..., Z_1, Z_0=1.0]
        Sk_array[i] = SofK(energy, T, L, i)
        dZdB_array[i] = Deriv_Partition((N - i) * energy, T, L)
        dSkdB_array[i] = Deriv_SofK(energy, T, L, i)

    Z_array[-1] = 1.0   # Naive implementation leaves these as "= L" otherwise
    Sk_array[0] = 1.0
    print(Z_array)
    print(Sk_array)
    print(dZdB_array)
    print(dSkdB_array)

    #denom = np.sum(Z_array)    # If these were bosons, this would work
    denom = 0.0
    numer = 0.0
    for i in range(Nplus):
        denom += eta**(i) * Sk_array[i] * Z_array[N - i]
        numer += eta**(i) * (dSkdB_array[i]*Z_array[N - i] + Sk_array[i]*dZdB_array[N - i])

    print(f"Energy using recursion: {numer / denom}")

    
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
