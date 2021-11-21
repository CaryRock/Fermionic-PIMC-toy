import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

with open("julia_beads_positions.csv") as file:
    matrix = list(csv.reader(file,delimiter=","))
file.close()

x = matrix[0]
y = np.linspace(1,20,20)

print(x)
print(y)

#plt.plot(x,y)
#plt.show()

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

plt.plot(x,y,'bo')
plt.show()
