import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt("dbg.txt")
plt.figure(figsize=(8, 5))
plt.plot(dat[:, 0], dat[:, 1], ".", ms=1)
plt.show()
