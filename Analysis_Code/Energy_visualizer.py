import numpy as np
import matplotlib.pyplot as plt
import sys
filename = sys.argv[1]

with open(filename) as f:
    z = f.read().splitlines()
f.close()

x = np.zeros(len(z))
y = np.zeros(len(z))

for i in range(len(z)):
    x[i] = i
    y[i] = float(z[i])

fig, ax = plt.subplots()

ax.set_title("System Energy")

ax.set_xlabel("MC step")
ax.set_ylabel("Energy")

ax.plot(x, y)
plt.show()
