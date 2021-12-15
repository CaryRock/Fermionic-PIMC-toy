import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    # Import the data file and read
    with open(sys.argv[1], "r") as f:
        data = np.loadtxt(f)
    f.close()

    # Get number of bins (len(data)), assume middle bin is 0 if odd #, middle two surround if even #
    N = len(data)
    MAX = max(data)

    # Subtract counts from each bin - for k = 0:N-1, bin[k] = data[N-1-k] - data[k]
    binned = np.zeros(int(N/2))
    for k in range(int(N/2)-1):
        binned[k] = data[N-1-k] - data[k]
        binned[k] /= MAX

    # Plot bin, should be ~= 0 for even, single spike in middle if odd
    fig = plt.figure()
    plt.plot(binned, 'b.')
    plt.title("Expectation values")
    plt.show()

if __name__ == "__main__":
    main()
