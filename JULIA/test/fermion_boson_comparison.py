import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    bosons = sys.argv[1]
    fermions = sys.argv[2]

    with open(bosons, "r") as bos:
        boson_data = np.loadtxt(bos)
    
    with open(fermions, "r") as fer:
        fermion_data = np.loadtxt(fer)
    
    bosons = boson_data[:,1]
    fermions = fermion_data[:,1]

    for i in range(len(bosons)):
        print(fermions[i] - bosons[i])
    
    temps = boson_data[:,0]#np.linspace(1, len(boson_data), len(boson_data))
    plt.plot(temps, bosons, 'b.', label="Boson <E>, N = 2")
    plt.plot(temps, fermions, 'r.', label="Fermion <E>, N = 2")
    plt.plot(temps, fermions - bosons, 'k.', label="Difference")
    plt.xlabel("T (K)")
    plt.ylabel("Energy (K)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
