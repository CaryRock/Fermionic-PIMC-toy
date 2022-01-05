#! /bin/usr/env python3
# Bins the error bars of a provided file. The file should have:
#       1.) The estimator desired
#       2.) The square of the estimator desired
#       3.) Sufficiently many steps such that at the desired level of binning, 
#           there will be >=30 data points being binned over
#           (i.e., binning in half 5 times is a reduction by a factor of 32)

import numpy as np
import sys
import docopt as docopt # Good to practice for future throw-away scripts


