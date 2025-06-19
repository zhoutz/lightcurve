import numpy as np
import matplotlib.pyplot as plt


dat1 = np.loadtxt("sd1.txt")

dat2 = np.loadtxt("apjlab5968/SD1a_test_IM.txt")
# dat2 = np.loadtxt("apjlab5968/SD1b_test_IM.txt")

plt.figure(figsize=(8, 5))

plt.plot(dat1[:, 0], dat1[:, 2], ".", ms=1)
# plt.plot(dat2[:, 0], dat2[:, 1], ".", ms=1)

plt.show()
