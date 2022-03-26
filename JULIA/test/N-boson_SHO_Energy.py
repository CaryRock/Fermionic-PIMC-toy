import numpy as np
import sys

def main():
    T_min = float(sys.argv[1])
    T_max = float(sys.argv[2])
    N = int(sys.argv[3])
    nPlus = N + 1
    temps = np.linspace(T_min, T_max, int(T_max - T_min) + 1)

    hbar = 1.0
    omega = 1.0
    delta = hbar * omega
    kB = 1.0
    
    energy = 0.0
    E_list = []
    
    print(temps)
    for t in temps:
        energy = 0.0
        for k in range(1, nPlus):
            energy += k / (np.exp(k * delta / t) - 1.0)
        energy += N/2.0

        E_list.append(energy)
        
        print(f"{t}\t{energy}")

if __name__ == "__main__":
    main()
