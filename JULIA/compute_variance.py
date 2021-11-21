import numpy as np
import sys
import matplotlib.pyplot as plt

def main():
    with open(sys.argv[1], "r") as f:
        X1_Data = f.read().splitlines()
    f.close

    with open(sys.argv[2], "r") as g:
        X2_Data = g.read().splitlines()
    g.close()

    N = len(X1_Data)
    u = np.zeros(N)
    v = np.zeros(N)
    w = np.zeros(N)

    for i in range(N):
        u[i] = float(X1_Data[i])
        v[i] = float(X2_Data[i])
        w[i] = v[i] - (u[i] * u[i])

    fig = plt.figure()
    plt.plot(w)
    plt.show()

if __name__ == "__main__":
    main()
