import sys
import numpy as np
import mpmath as mpm
import matplotlib.pyplot as plt
from math import exp

def main():
    # temp_dat.txt
    N = 0
    numArgs = len(sys.argv)
    print(f"sys.argv = {sys.argv}")
    with open(sys.argv[1],"r") as f:
        temp = f.read().splitlines()
        N = len(temp)
    f.close()

    # x1_dat.txt
    with open(sys.argv[2],"r") as g:
        x1_raw = g.read().splitlines()
    g.close()

    # x2_dat.txt
    with open(sys.argv[3],"r") as h:
        x2_raw = h.read().splitlines()
    h.close()
    
    # x2_err.txt
    if(numArgs == 5):
        with open(sys.argv[4],"r") as i:
            x2_err = i.read().splitlines()
        i.close()

    N = len(temp)
    
    T       = np.zeros(N)
    x1      = np.zeros(N)
    x2      = np.zeros(N)
    for j in range(N):
        T[j] = float(temp[j])
        x1[j] = float(x1_raw[j])
        x2[j] = float(x2_raw[j])
    
    if(numArgs == 5):
        x2err = np.zeros(N)
        for k in range(N):
            x2err[j] = float(x2_err[j])

    x   = np.linspace(-2.0,2.0,N)
    xp  = np.linspace(-2.0,2.0,N)

    rho = np.zeros(N)

    l02 = 1.0               # hbar/(m omega)
    for i in range(N):
        eps = 1.0 / T[i]           # beta*hbar*omega
        sinhEps = mpm.sinh(eps)
        cothEps = mpm.coth(eps)
        denom = 1/np.sqrt(mpm.coth(eps/2)*np.pi*l02)
#        rho[i] = denom * exp( 1/(2*l02)* (-cothEps * (x1[i]*x1[i] + x1[j]*x1[j]) + 2*x[i]*x1[j]/sinhEps))
        rho[i] = 1/2 * mpm.coth(1/(2*T[i]))

    print(f"T = {T}\nx1 = {x1}\nx2 = {x2}")
    fig = plt.figure()
    plt.plot(T,rho,'ro', label="rho")
    if(numArgs == 5):
        plt.errorbar(T,x2,x2err,fmt='bo', label="x2")
    else:
        plt.plot(T,x2,'bo',label="x2")
    plt.xlabel("Temp (K)")
    plt.ylabel("E (K)")
    plt.legend()
    plt.yscale('log')
    plt.show()

if __name__ == "__main__":
    main()
