import numpy as np
import matplotlib.pyplot as plt

size = [4, 12, 60, 660, 1980, 2448, 3960, 5544, 12240]
time = [1.9e-05, 4.1e-05, 0.000621, 0.260119, 9.47873, 16.9626 ,66.7833, 184.291, 2056.3]

plt.loglog(size, time, 'o-')
plt.xlabel("Number of states")
plt.ylabel("Time [s]")
plt.show()
