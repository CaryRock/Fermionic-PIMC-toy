#! /usr/bin/env python3
# TODO: (maybe?) IMPLEMENT A VERSION THAT IS COMPATIBLE WITH PYTHON 2.7
# This script reads in a .dat file for ce-lineardensity-* and outputs the
# difference between the input data and the theoretical curve it should match.

import sys
import numpy as np
import scipy.integrate
import mpmath as mpm
from math import exp
import matplotlib.pyplot as plt
import argparse
import pathlib
import glob

def create_parser():
    parser = argparse.ArgumentParser(description='First iteration of a script \
to plot the difference between the theoretical and Monte Carlo-generated data \
for a single particle in a simple harmonic oscillator in 1-D.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-T", "--temp", type=float, help="Required. Temperature of data - single \
data input, so single temperature value (in K).", required=True)
    
    parser.add_argument("--data", type=pathlib.Path, help="Required. File \
containing X1 data", required=True)

    parser.add_argument("--figName", type=str, default='', help="Optional. Name under which to save output figure.")

    # TODO: IMPLEMENT FUNCTIONALITY TO PARSE ALL INPUT FILES AT SAME TIME AND
    # OUTPUT A SINGLE GRID PLOT OF THEM ALL (ALL AT ONCE)
    # That would be glob.glob("ce-{estimator}-...")

    return parser

def main(argv=None):
### Main Variabls, Constants, etc. ############################################
    if argv is None:
        argv = sys.argv

    parser = create_parser()
    args = parser.parse_args(argv[1:])

    T           = args.temp
    alpha       = 1.0   # hbar * omega / k_B
    lam         = 1/2
    lO2         = lam / np.sqrt(alpha/2)
    eps         = np.sqrt( alpha/2 ) / T * 1.0 # hbar/T * sqrt(alpha/2)
    
    sinhEps     = mpm.sinh(eps)
    cothEps     = mpm.coth(eps)
    denom       = 1/np.sqrt(mpm.coth(eps/2) * np.pi / lO2)
    
    p = [denom, 1/(2 * lO2), cothEps, sinhEps]
    figName     = args.figName
    print(figName)
### Read the data file ###
    with open(args.data) as f:
        lines = (line for line in f if not line.startswith('#'))
        data = np.loadtxt(lines)  # This is a stupid name - change that TODO
    f.close()
    
    n = len(data)

    print(f"The number of data lines in the file is: {n}")

### Create variables and solve for the exact answer; prepare the plot #########
    x_a = data[0][0]  # These should be the values of the first and last rows of the input data file
    x_b = data[n-1][0]
    
    delta = data[1][0] - data[0][0]
   
    counts = np.zeros(n)            # The data curve
    u   = np.linspace(x_a, x_b, n)  # Setup the x-axis values
    for i in range(n):
        u[i] += delta/2             # To align to center of bins
        counts[i] = data[i][1]        # Could do without this, but separats out count values nicely
    v   = np.copy(u)
    rho = np.zeros(n)               # Theoretical curve

    for i in range(n):
        rho[i] = p[0] * exp( p[1] * (-p[2] * (u[i]*u[i] + v[i]*v[i]) + 2*u[i]*v[i]/p[3]))
        # Compute that curve

    # Normalize the theoretical curve
    norm = scipy.integrate.simps(rho,u)
    rho /= norm
    
    # Normalize the data curve
    norm = scipy.integrate.simps(counts,u)
    counts /= norm

    diff = np.zeros(n)
    for j in range(n):
        diff[j] = rho[j] - counts[j]
    
    chi2 = 0.0
    count = 0
    count2 = 0
    for k in range(n):
        chi2 += diff[k] * diff[k]
        count2 += 1

    chi2 /= n

    print(f"chi^2 = {chi2}")

    plt.figure()
    plt.axhline(y=0,alpha=0.5)
    plt.plot(u, rho, 'r--', alpha=0.5, label="Theoretical curve")
    plt.plot(u, counts, 'k--', alpha=0.5, label="MC data")
    plt.plot(u, diff, 'b.', label="Difference between Theory and MC data")
    plt.suptitle(f"Lineardensity @ T = {T} K")
    plt.title(f"File={args.data}")
    plt.xlabel("Position")
    plt.ylabel("Arb.")
    plt.figtext(0.2, 0.85, r"$\chi^2 = {}$".format(chi2))
    plt.legend()
    
    if args.figName != '':
        plt.savefig(f"./{figName}.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()

