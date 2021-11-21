import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def autocorrelation_function(time, scale, autocorrelation_time):
    return scale*np.exp(-time/autocorrelation_time)

def autocorrelation(data):
    N = data.shape[0]
    _autocorrelation = np.zeros(N)
    for dt in range(N-1):
        denom = np.mean(data[:N - dt]**2) - np.mean(data[:N - dt])**2
        num = np.mean(data[:N - dt]*data[dt:]) - np.mean(data[:N - dt]) * \
            np.mean(data[dt:])
        _autocorrelation[dt] = num/denom
    return _autocorrelation

def main():
    filename = sys.argv[1]
    equil = 0
    observ = 0
    fTemp = 0

    # Import data to be analyzed
    with open(filename, "r") as f:
        dat = f.read().splitlines()
    
    equil = int(dat[0])
    observ = int(dat[1])
    fTemp = float(dat[2])

    x = np.zeros(len(dat)-3)
    N = len(x)
    for i in range(N):
        x[i] = float(dat[i+3])

    # Instantiate the relevant variables - _autocorrelation, etc.
    _autocorrelation = autocorrelation(x)
    time_separation = np.arange(_autocorrelation.shape[0])

    # Use the magic of curve_fit to determine the correlation length
    popt, perr = curve_fit(autocorrelation_function, time_separation, _autocorrelation)

    # Results
    print("t_auto: {}".format(popt[1]))
    print("5t_auto: {}".format(5*popt[1]))

    # Plot
    fig, ax = plt.subplots(figsize=(8,4.5), dpi=120, constrained_layout=True)
     
    ax.plot(time_separation, _autocorrelation)
    ax.plot(time_separation, autocorrelation_function(time_separation, *popt))
    #ax.set_xlim((0,100))    
    plt.show()

if __name__ == "__main__":
    main()
