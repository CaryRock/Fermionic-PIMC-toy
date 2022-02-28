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

    temps = np.linspace(1, len(boson_data), len(boson_data))
    plt.plot(temps, boson_data, 'b.', label="Boson <E>, N = 2")
    plt.plot(temps, fermion_data, 'r.', label="Fermion <E>, N = 2")
    plt.plot(temps, fermion_data - boson_data, 'k.', label="Difference")
    plt.xlabel("T (K)")
    plt.ylabel("Energy (K)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
