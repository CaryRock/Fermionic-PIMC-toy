import sys
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mpm

def main():
    with open(sys.argv[1], "r") as f:
        X2_temp = f.read().splitlines()
    f.close()

    with open(sys.argv[2], "r") as g:
        X2_data = g.read().splitlines()
    g.close()

    with open(sys.argv[3], "r") as h:
        X2_err = h.read().splitlines()
    h.close()

    N = len(X2_temp)
    x = np.zeros(N)
    y = np.zeros(N)
    err = np.zeros(N)

    for i in range(N):
        x[i] = float(X2_temp[i])
        y[i] = float(X2_data[i])
        err[i] = float(X2_err[i])
    
    z = np.polyfit(x,y,1)
    p = np.poly1d(z)
    
    a = np.linspace(0.5,1.5,101)
    b = np.linspace(0.5,1.5,101)
    for i in range(100):
        a[i] = 1/2*mpm.coth(1/(2*b[i]))

    fig = plt.figure()
#    plt.plot(x, y, 'bo')
    plt.plot(b,a, 'go')
    plt.plot(x,p(x), "r--",label="0.8363x + 0.3361")
    plt.errorbar(x,y,err,fmt='bo')
    print(f"y={z[0]}x + {z[1]}")
    plt.title(r"$x^2$ vs Temp; N = 250 runs; N_MC = 100,000")
    plt.xlabel('Temp (K)')
    plt.ylabel(r'$<X^2>$')
#    plt.legend(p,"0.8363x + 0.3361")
    plt.show()



if __name__ == "__main__":
    main()
