import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def line_test(x, m, b):
    return m*x + b

def main():
    print(f"sys.argv: {sys.argv}")
    print(f"Length of sys.argv: {len(sys.argv)}")
    print(f"sys.argv[2]: {sys.argv[2]}")

    zeros = []
    for i in range(len(sys.argv)):
        if i == 0: continue
        with open(sys.argv[i], "r") as f:
            data = f.read().splitlines()
        f.close()
        
        sum = 0.0
        for j in range(len(data)):
            if j == 0: continue
            sum += float(data[j])

        zeros.append(sum)

        #print(f"{zeros}")

    fig = plt.figure()
    ax = plt.gca()
    xLabels = np.arange(0, 2200, 200)
    xLabels[0] = 100
    xticks = np.arange(0, len(zeros), 1)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xLabels)
    x_dat = np.arange(0, len(zeros), 1)
    
    params, params_covar = optimize.curve_fit(line_test, x_dat, zeros, p0=[0.01, 0])

    plt.plot(zeros,'b.', label="Expectation Value Data")
    plt.plot(x_dat, line_test(x_dat, params[0], params[1]), 'r--',label="Fitted line")
    plt.ylabel("< x > (L)")
    plt.xlabel("# of Beads")
    plt.legend()
#    plt.suptitle("")
    plt.title("< x > vs # of Beads in that Trial")
    plt.show()
    
if __name__ == "__main__":
    main()
