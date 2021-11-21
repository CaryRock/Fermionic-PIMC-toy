import sys
import numpy as np
import matplotlib.pyplot as plt
import mpmath as mpm

def Exact(temp):
    return 1/2*mpm.coth(1/(2*temp))

def Average(x,N):
    tot = 0
    for i in range(N):
        tot += x[i]

    return tot/N

def main():
    filename = sys.argv[1]
    equil = 0
    observ = 0
    fTemp = 0

    with open(filename,"r") as f:
        dat = f.read().splitlines()
    f.close()

    equil = int(dat[0])
    observ = int(dat[1])
    fTemp = float(dat[2])

    x = np.zeros(len(dat)-3)
    N = len(x)

    for i in range(N):
        x[i] = float(dat[i+3])

    print(f"Exact energy: {Exact(fTemp)}")
    print(f"Average energy: {Average(x,N)}")

if __name__ == "__main__":
    main()
