import sys
import numpy as np
import scipy.integrate
import mpmath as mpm
from math import exp
import matplotlib.pyplot as plt
import argparse
import pathlib
import re
import concurrent.futures
import multiprocessing as mlt   # Gratuitous to have both; am nub

def create_parser():
    parser = argparse.ArgumentParser(description='Temperature-dependent \
            graphic of x1 and/or x2 data.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("temp", type=float, help="Temperature of the data")
    parser.add_argument("-x1", "--x1Data", type=pathlib.Path, help="File \
            containing x1 data.")
    parser.add_argument("-x2", "--x2Data", type=pathlib.Path, help="File \
            containing x2 data.")
    parser.add_argument("-b1", "--binnedX1", type=pathlib.Path, help="File \
            containing binned X1 data.")
    parser.add_argument("-b2", "--binnedX2", type=pathlib.Path, help="File \
            containing binned X2 data.")
    parser.add_argument("-prod", "--production", type=pathlib.Path, help="File \
            containing production PIMC data")
    parser.add_argument("-nbins", "--numBins", type=int, default=51, 
        help="Number of bins to use in analysis.")
    parser.add_argument("-lim", type=float, default=3, help="min & max values \
            of exact curve.")
    parser.add_argument("-err", "--errorBars", action='store_true', help="Add \
            error bars on the read data.")
    parser.add_argument("--binRho", action="store_true", help="Bin the output \
of rho to compare to a binned input file.")
    parser.add_argument("--title", type=pathlib.Path, help="Because I'm lazy \
            and this is an easy way to get the title of the \
            graphs consistently.")
    return parser

def make_plot(data, numBins, errorBars, legLab=""):
    M = len(data)
    s = np.zeros(M)
    if type(data[0]) is not float:
        for a in range(M):
            s[a] = float(data[a])
            # If data's already float, then this makes no difference

    x_a = min(s)
    x_b = max(s)

    binDelta = (x_b - x_a)/numBins
    x = np.linspace(x_a, x_b, numBins)
    y = np.zeros_like(x)

    binLoc = 0
    for b in range(M):
        binLoc = int(s[b] / binDelta + x_b)
        if binLoc >= numBins: binLoc = numBins - 1
        elif binLoc < 0: binLoc = 0
        y[binLoc] += 1

    norm = scipy.integrate.simps(y,x)
    y /= norm
    print(f"Old integral value of data: I = {norm}")
    new = scipy.integrate.simps(y,x)
    print(f"Normalized value of integral: I = {new}")
    yStdErr = np.zeros_like(y)
    for c in range(numBins):
        yStdErr[c] = np.sqrt(y[c]/norm)

    sum = 0
    for d in range(numBins):
        sum += y[d]

    sum /= binDelta
    
    if errorBars:
        plt.errorbar(x,y, yerr=yStdErr, fmt='b.', label=legLab)
    else:
        plt.plot(x, y, 'b.', label=legLab)

def make_binned_plot(data, numBins, errorBars, rhoMax):
    M = len(data)
    s = np.zeros(M)
    if type(data[0]) is not float:
        for a in range(M):
            s[a] = float(data[a])

    x_a = -3.0
    x_b = 3.0
    binDelta = (x_b - x_a)/numBins

    x = np.linspace(x_a, x_b, M)
    norm = scipy.integrate.simps(s, x)
    s /= norm
    
    print(f"Old integral value of data: I = {norm}")
    new = scipy.integrate.simps(s,x)
    print(f"Normalized value of integral: I = {new}")
    yStdErr = np.zeros_like(s)
    for c in range(numBins):
        yStdErr[c] = np.sqrt(s[c]/norm)

    if errorBars:
        plt.errorbar(x,s, yerr=yStdErr, fmt='b.')
    else:
        plt.plot(x, s, 'b.', label="Norm'd, binned")
    
    sPrime = np.copy(s)
    maxNormS = max(s)
    print(f"max(s) = {maxNormS}\nmax(rho) = {rhoMax}\n")
    for e in range(len(s)):
        sPrime[e] *= rhoMax / maxNormS
    
    #norm = scipy.integrate.simps(sPrime,x)
    #sPrime /= norm
    plt.plot(x, sPrime, 'k.',label="rescaled binned") 


# TODO: These two are the same function. Reduce, reuse, recycle.
# def read_Dat(filename, numBins, errorBars, binned, legLab):
def x1_data(filename, numBins, errorBars, binned, legLab, rhoMax):
    with open(filename,"r") as f:
        raw_x1_data = np.loadtxt(f)
    f.close()

    if not binned:
        make_plot(raw_x1_data, numBins, errorBars, legLab)
    else:
        make_binned_plot(raw_x1_data, numBins, errorBars, rhoMax)

def x2_data(filename, numBins, errorBars, binned):
    with open(filename,"r") as g:
        raw_x2_data = g.read().splitlines()
    g.close()

    if not binned:
        make_plot(raw_x2_data, numBins, errorBars, legLab)
    else:
        make_binned_plot(raw_x2_data, numBins, errorBars, rhoMax)

def parse_prod(name):
    with open(name, "r") as f:
        raw_read = np.loadtxt(f)
    f.close()

    N = len(raw_read)
    x = raw_read[2:,0]
    y = raw_read[2:,1]
    yError = raw_read[2:,2]

    return x, y, yError

    
def Rho(ui, v, p, n):
    rhoJ = np.zeros(n)
    for j in range(n):
        rhoJ[j] = p[0] * exp( p[1] * (-p[2] * (ui*ui + v[j]*v[j]) + 2*ui*v[j]/p[3]))
    return rhoJ

def main(argv=None):
### Main Variables, Constants, etc. ###########################################
    if argv is None:
        argv = sys.argv

    parser = create_parser()
    args = parser.parse_args(argv[1:])

    T           = args.temp
    alpha       = 1.0   # hbar * omega / k_B
    lam         = 1/2
    lO2         = lam / np.sqrt(alpha/2)
    eps         = np.sqrt( alpha/2 ) / T * 1.0  # hbar/T * sqrt(alpha/2)

    sinhEps     = mpm.sinh(eps)
    cothEps     = mpm.coth(eps)
    denom       = 1/np.sqrt(mpm.coth(eps/2) * np.pi / lO2)
    
    p = [denom, 1/(2*lO2), cothEps, sinhEps]
    numBins     = args.numBins


### Create variables and solve for the exact answer; prepare plot #############
    x_a = -args.lim
    x_b = args.lim
    if args.binRho == True:
        n = 100000
    else:
        n = 100

    u   = np.linspace(x_a, x_b, n)  # Basically x in exact solution
    v   = np.linspace(x_a, x_b, n)  # Basically x' in exact solution
    rho = np.zeros([n,n])

    #for i in range(n):
    #    rho[i] = Rho(u[i], v, p, n)
#    with concurrent.futures.ProcessPoolExecutor() as executor:
#        for i in range(n):
#            rho[i] = Rho(u[i], v, p, n)
#
#    for j in range(n):  # See about using multiprocess to speed this up
#        rho[j] = rho[j].result()
    
#        for j in range(n):
#            rho[i][j] = denom * exp( 1/(2*lO2) * (-cothEps * (u[i]*u[i] + v[j]*v[j]) + 2*u[i]*v[j]/sinhEps))
#            rho[i][j] = p[0] * exp( p[1] * (-p[2] * (u[i]*u[i] + v[j]*v[j]) + 2*u[i]*v[j]/p[3]))
    for i in range(n):
        rho[i][i] = p[0] * exp( p[1] * (-p[2] * (u[i]*u[i] + v[i]*v[i]) + 2*u[i]*v[i]/p[3]))
    
    traceSum = 0.0
    rho2 = 0.0
    rhoMax = 0.0
    if args.binRho == False:
        trace = np.zeros_like(u)
        for l in range(n):
            trace[l] = rho[l][l]   #(x_b - x_a)*rho[l][l]/n
    
        rhoNorm = scipy.integrate.simps(trace,u)
    
        print(f"rhoNorm = {rhoNorm}")
        trace /= rhoNorm
        rhoMax = max(trace)
        print(f"Normal'd value of theo. integral: I = {scipy.integrate.simps(trace,u)}")
        fig = plt.figure()
        plt.plot(u, trace, 'r--', label="theory")
    else:
        _wid = 100
        L = int(n/_wid)
        w = np.linspace(x_a, x_b, L)
        trace = np.zeros_like(w)
        for l in range(L):
            for m in range(_wid):
                trace[l] += rho[l*_wid+m][l*_wid+m]

        rhoNorm = scipy.integrate.simps(trace,w)
        rhoMax = max(trace)
        print(f"rhoNorm = {rhoNorm}")
        trace /= rhoNorm
        print(f"Normal'd value of theo. integral: I = {scipy.integrate.simps(trace,w)}")
        fig = plt.figure()
        plt.plot(w, trace, 'r--', label="Binned theory")


### X1 data-dependent functions ###############################################
    do_x1 = False
    if args.x1Data != None or args.binnedX1 != None:
        do_x1 = True
        if not args.binnedX1: 
            x1_data(args.x1Data, numBins, args.errorBars, False, legLab, 1)
        else:
            x1_data(args.binnedX1, numBins, args.errorBars, True, "test", rhoMax)


### X2 data-dependent functions ###############################################
    do_x2 = False
    if args.x2Data != None or args.binnedX2 != None:
        do_x2 = True
        if not args.binnedX2:
            x2_data(args.x2Data, numBins, args.errorBars, False, legLab, 1)
        else:
            x2_data(args.x2Data, numBins, args.errorBars, True, "testX2", rnoMax)


### Account for data being from Production_Pimc code ##########################
    do_prod = False
    if args.production != None:
        do_prod = True
        x, y, yError = parse_prod(args.production)
#        make_plot(data, numBins, args.errorBars)
        delt = x[1] - x[0]
#        plt.plot(x, y, 'k.')
        plt.bar(x, y, yerr=yError, width=delt, edgecolor='k')


### Error on no input and plot results ########################################
    if (do_x1 == False) and (do_x2 == False) and (do_prod == False):
        print("Please run with data for x1 and/or x2 and/or production data\n")
        exit()

### Plot title and file-name ##################################################
    plt.suptitle(fr"Counts and $\rho$ vs Position at {T} K")
    if args.title is None:
        plt.title(f"file={args.x1Data}", y=1)
    else:
        plt.title(f"file={args.title}", y=1)
    plt.xlabel("x (L)")
    plt.ylabel(r"Counts (arb.) or $\rho$ (1/L)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
