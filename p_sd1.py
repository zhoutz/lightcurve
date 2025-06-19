import numpy as np
import matplotlib.pyplot as plt


dat1 = np.loadtxt("sd1.txt")

# dat2 = np.loadtxt("apjlab5968/SD1a_test_IM.txt")
# dat2 = np.loadtxt("apjlab5968/SD1b_test_IM.txt")
dat2 = np.loadtxt("apjlab5968/SD1c_test_IM.txt")

k = np.gradient(dat1[:, 0], dat1[:, 0] + dat1[:, 2])

plt.figure(figsize=(8, 5))

# plt.plot(dat1[:, 0], dat1[:, 1], ".", ms=1)
plt.plot(dat1[:, 0] + dat1[:, 2], dat1[:, 1] * k, ".", ms=1)
plt.plot(dat2[:, 0], dat2[:, 1], ".", ms=1)

plt.show()
