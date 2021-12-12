#! /usr/bin/env python3
# Test to return only the desired estimator files

import numpy as np
import matplotlib.pyplot
import sys
import os
import glob

def get_binned_error(mc_data):
    '''Get the standard error in mc_data and return neighbor averaged data.'''
    N_bins = mc_data.size
    delta = np.std(mc_data)/np.sqrt(N_bins)

    start_bin = N_bins % 2
    binned_mc_data = 0.5*(mc_datapstart_bin::2] + mc_data[start_bin + 1::2])
    
    #This reduces the length of binned_mc_data by ~2x (-> 1/2 original)
    return delta, binned_mc_data

def main():
#    path =     # Glob by default looks in current directory
    files = glob.glob("ce-estimator-*.dat")

    for i in range(len(files)):
        print(files[i])
    
    # M is from this project - entirely different data. Obviously should be
    # rewritten for the current affiar
    nLevs = in(np.log2(M[skip:,0].size/4))+1
if __name__ == "__main__":
    main()
