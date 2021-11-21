import sys
import numpy as np
import mpmath as mpm
import matplotlib.pyplot as plt
from math import exp

def main():
    with open(sys.argv[1],"r") as f:
        t = len(f.readlines())
        for i in range(t):
            temp = f.read().splitlines()
    
    for j in range(len(temp)):
        T = float(temp[j])

    l02 = 1.0               # hbar/(m omega)

    x   = np.linspace(-2.0,2.0,100)
    xp  = np.linspace(-2.0,2.0,100)
    N = len(x)

    rho = np.zeros([N,N])

    for i in range(N):
        eps = 1.0 / T           # beta*hbar*omega
        sinhEps = mpm.sinh(eps)
        cothEps = mpm.coth(eps)
        denom = 1/np.sqrt(mpm.coth(eps/2)*np.pi*l02)
        rho[i] = denom * exp( 1/(2*l02)* (-cothEps * (x[i]*x[i] + xp[j]*xp[j]) + 2*x[i]*xp[j]/sinhEps))
    
    fig = plt.figure()
    plt.plot(x,rho,'r--')
    plt.xlabel("x")
    plt.ylabel("rho")
    plt.show()

if __name__ == "__main__":
    main()
