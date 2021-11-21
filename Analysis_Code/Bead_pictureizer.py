import sys                      # Command line input
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd             # Historgram

filename = sys.argv[1] 

with open(filename) as f:
    z = f.read().splitlines()

x = np.zeros(len(z))

for i in range(len(z)):
    x[i] = 100*float(z[i])

y = np.zeros(len(x))
for i in range(len(x)):
    y[i] = i
'''
fig, ax = plt.subplots()

ax.set_title("Bead 0 Time Slice positions")

ax.set_xlabel("position")
ax.set_ylabel("time slice #")

xMax = np.floor(max(x)) + 1
ax.set_xlim(-xMax, xMax)
ax.set_ylim(0,40)

ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))

ax.plot(x, y,'bo')
plt.show()

resolution = 100
numBins = 2*int(xMax)*resolution

print(len(z))
slices = np.linspace(0, len(z), numBins+1, dtype=int)
binnedPositions = np.add.reduceat(x, slices[:-1])/np.diff(slices)

fig, ax = plt.subplots()
ax.set_title("Bead 0 Time Slice positions - histogram")
ax.set_xlabel("Bins")
ax.set_ylabel("# of particles")

ax.plot(binnedPositions,'ro')
plt.show()
'''

plt.hist(x,bins=range(-300,300,5))
plt.show()
