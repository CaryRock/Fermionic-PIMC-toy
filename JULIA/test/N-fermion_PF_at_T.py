import numpy as np
import sys

"""
Consider: what are the high and low temperature limits that should be considered?
In the high-temperature limit, the energies should basically approach those of
the classical, distinguishable particles (governed by Boltzmann statistics). In
the low-temperature limit, they should definitely be quantum, different. The job
of this code: to connect the dots.
"""

def Single_Energy(j):
    return j + 0.5

class PF:
    def __init__(self, T, L, energy):
        self.T = T                  # The temperature
        self.L = L                  # The cutoff length
        self.Z = 0.0                # The value of the partition function
        self.dZdB = 0.0             # The value of the derivative of Z wrt beta
        self.energy = energy        # The energy spectrum

    def Partition(self, zero=0):
        if (zero != 0):
            sum = 0.0
    
            for i in range(self.L):
                sum += np.exp(-self.energy[i] / self.T)
    
            self.Z = sum

        else:
            self.Z = 1.0

    def Deriv_Partition(self, zero=0):
        if (zero != 0):
            sum = 0.0

            for i in range(self.L):
                sum += self.energy[i] * np.exp(-self.energy[i] / self.T)

                self.dZdB = sum

        else:
            self.dZdB = np.sum(self.energy)

class Sk:
    def __init__(self, T, L, energy, multiple):    # "multiple" should never = 0
        self.T = T                  # The temperature
        self.L = L                  # The cutoff length
        self.k = multiple           # The "k" of S(k)
        self.Sk = 0.0               # The value of the S(k) function
        self.dSkdB = 0.0            # The derivative of S(k) wrt beta
        self.energy = energy        # The energy spectrum previously computed
    
    def SofK(self):
        sum = 0.0
    
        for i in range(self.L):
            sum += np.exp(-self.energy[i] * self.k / self.T)
    
        return sum
    
    def Deriv_SofK(self):
        sum = 0.0
    
        for i in range(self.L):
            sum += self.energy[i] * np.exp(-self.energy[i] * self.k / self.T)
    
        return (-self.k * sum)

def main():
    T = float(sys.argv[1])
    N = int(sys.argv[2])
    Nplus = N + 1

    hbar = 1.0
    w = 1.0
    kB = 1.0
    L = int(10.0 * T)  # Cutoff - the number of energy levels to compute
    # When is there a sufficient number of energy levels taken? When E/T >> 1 => e^-(E/T) approx 1
    eta = -1
    energy = np.zeros(L)
    for i in range(L):
        energy[i] = Single_Energy(i)
    
    """
    Since the energy is now variable depending on the temperature that the PF
    is being evaluated at, should bring back the class-based definition of Z 
    and have each Z record its own energy spectrum, Sk_array, and the corresponding
    derivatives.
    """
    '''
    Z_array = [PF(T, L, energy) for i in range(Nplus)] # Array to store PFs - N+1 because 0 -> N
    dZdB_array = [1.0 for i in range(Nplus)]
    Sk_array = [Sk(T, L, energy, i) for i in range(1, Nplus)]
    dSkdB_array = [1.0 for i in range(Nplus)]

    for i in range(Nplus):
        Z_array[i] = Partition((N - i) * energy, T, L)   # Produces Zarray[Z_N, Z_N-1, ..., Z_1, Z_0=1.0]
        Sk_array[i] = SofK(energy, T, L, i)
        dZdB_array[i] = Deriv_Partition((N - i) * energy, T, L)
        dSkdB_array[i] = Deriv_SofK(energy, T, L, i)

    Z_array[-1].Partition(0)
    #Sk_array[0].SofK(0) = 1.0  # There shouldn't be a "S(0)" term
    #print(Z_array)
    #print(Sk_array)
    #print(dZdB_array)
    #print(dSkdB_array)

    #denom = np.sum(Z_array)    # If these were bosons, this would work
    denom = 0.0
    numer = 0.0
    for i in range(1,Nplus):
        denom += eta**(i) * Sk_array[i] * Z_array[N - i]
        numer += eta**(i) * (dSkdB_array[i]*Z_array[N - i] + Sk_array[i]*dZdB_array[N - i])

    print(f"Energy using recursion: {numer / denom}")
    '''
    
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
    
    denom = sum(Z1) - sum(Z2)
    numer = 0.0
    for i in range(L):
        for j in range(L):
            numer += (energy[i] + energy[j]) * np.exp(-(energy[i] + energy[j])/T)

        numer -= 2.0 * energy[i] * np.exp(-2 * energy[i] / T)


    E = numer / denom

    print(f"T = {T}\tE = {E}")

if __name__ == "__main__":
    main()
