import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate
import mpmath as mpm
from math import exp

def main():
    T       = float(sys.argv[1])
    alpha   = 1.0
    lam     = 1/2
    lO2     = lam/np.sqrt(alpha/2)
    eps     = np.sqrt(alpha/2) / T

    sinhEps = mpm.sinh(eps)
    cothEps = mpm.coth(eps)
    denom   = 1/np.sqrt(mpm.coth(eps/2) * np.pi / lO2)
    p = [denom, 1/(2*lO2), cothEps, sinhEps]

    with open(sys.argv[2], "r") as f:
        next(f)
        next(f)
        data1 = np.loadtxt(f)
    f.close

    with open(sys.argv[3], "r") as g:
        next(g)
        next(g)
        data2 = np.loadtxt(g)
    g.close()
    
    x = data1[:,0]
    yPrime = data1[:,1]
    y = yPrime / scipy.integrate.simps(yPrime, x)
    yError = y / np.sqrt(len(y))
    u = data2[:,0]
    v = data2[:,1]
    vError = data2[:,2]
    
    print(f"Max of y: {max(y)}")
    print(f"Index of max: {np.argmax(max(y))}")
    print(f"Location of maximum of binned data: [{x[np.argmax(max(y))]}, {max(y)}]")

    x_a = min(x)
    x_b = max(x)
    n = 1000
    a = np.linspace(x_a, x_b, n)
    b = np.linspace(x_a, x_b, n)
    rho = np.zeros(n)
    
    for i in range(n):
        rho[i] = p[0] * exp( p[1] * (-p[2] * (a[i]*a[i] + b[i]*b[i]) + 2*a[i]*b[i]/p[3] ))
    rhoNormed = rho/ scipy.integrate.simps(rho, a)

### Plotting Stuff ############################################################
    fig = plt.figure()
    plt.plot(a, rhoNormed, 'k--', label="Theory")
    plt.errorbar(x,y, yerr=yError, fmt='b.', label="Julia")
    plt.errorbar(u,v, yerr=vError, fmt='r.', label="Production")
    plt.legend()
    plt.suptitle("Counts vs Position")
    plt.title("Prod: Benchmark     Julia: -T 0.5 -N 1 -J 200 -E 10,000 -O 5000 -S 1000")
    plt.xlabel("X (A)")
    plt.ylabel("Counts (arb.)")
    plt.show()

if __name__ == "__main__":
    main()
