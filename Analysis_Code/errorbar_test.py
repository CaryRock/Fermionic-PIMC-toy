import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
x = np.arange(10)
y = 2.5 * np.sin(x / 20 * np.pi)
yerr = np.linspace(0.05,0.2,10)

plt.errorbar(x, y, yerr=yerr, label="Both limits (the default)")
plt.hlines(3,x[0],x[9])
plt.legend(loc="lower right")

fig.savefig("test.png")
plt.show()
