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
import re

def create_parser():
    parser = argparse.ArgumentParser(description='First iteration of a script \
to plot the difference between the theoretical and Monte Carlo-generated data \
for a single particle in a simple harmonic oscillator in 1-D.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-T", "--temp", type=float, help="Required. Temperature of data - single \
data input, so single temperature value (in K).", required=True)
    parser.add_argument("--data", type=pathlib.Path, help="Required. File \
containing X1 data", required=True)


    # TODO: IMPLEMENT FUNCTIONALITY TO PARSE ALL INPUT FILES AT SAME TIME AND
    # OUTPUT A SINGLE GRID PLOT OF THEM ALL (ALL AT ONCE)

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
### Read the data file ###
    with open(args.data) as f:
        lines = (line for line in f if not line.startswith('#'))
        FH = np.loadtxt(lines)  # This is a stupid name - change that TODO
    f.close()

    print(f"The number of data lines in the file is: {len(FH)}")

    n = len(FH)
    #print(f"x_a = {FH[0][0]}")
    #print(f"x_b = {FH[n-1][0]}")

### Create variables and solve for the exact answer; prepare the plot #########
    x_a = FH[0][0]  # These should be the values of the first and last rows of the input data file
    x_b = FH[n-1][0]

    delta = FH[1][0] - FH[0][0]
   
    counts = np.zeros(n)
    u   = np.linspace(x_a, x_b, n) #Currently, has 200 data lines and 2 comment lines
    for i in range(n):
        u[i] += delta/2
        counts[i] = FH[i][1]
    v   = np.copy(u)
    rho = np.zeros(n)

    for i in range(n):
        rho[i] = p[0] * exp( p[1] * (-p[2] * (u[i]*u[i] + v[i]*v[i]) + 2*u[i]*v[i]/p[3]))

    norm = scipy.integrate.simps(rho,u)
    rho /= norm
    norm = scipy.integrate.simps(counts,u)
    counts /= norm

    diff = np.zeros(n)
    for j in range(n):
        diff[j] = rho[j] - counts[j]
    print(f"rho[100] = {rho[100]}")
    print(f"counts[100] = {counts[100]}")
    print(f"diff[100] = {diff[100]}")

    

    fig = plt.figure()
    plt.axhline(y=0)
#    plt.plot(u, rho, 'r--', label="Theoretical curve")
    plt.plot(u, diff, 'b.', label="Difference between Theory and MC data")
    plt.title(f"T = {T}    File={args.data}")
    plt.legend()
    plt.show()

    

if __name__ == "__main__":
    main()

