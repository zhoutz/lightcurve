import numpy as np
import matplotlib.pyplot as plt
import sys


s = sys.argv[1]

dat1 = np.loadtxt("sd1.txt")
dat2 = np.loadtxt(f"apjlab5968/SD1{s}_test_IM.txt")

plt.figure(figsize=(8, 5))

# plt.subplot(2,1,1)

if True:
    plt.plot(dat1[:, 0], dat1[:, 1], ".", ms=1)
else:
    k = np.gradient(dat1[:, 0], dat1[:, 0] + dat1[:, 2])
    plt.plot(dat1[:, 0] + dat1[:, 2], dat1[:, 1] * k, ".", ms=1)

plt.plot(dat2[:, 0], dat2[:, 1], ".", ms=1)

if True:
    yhat = np.interp(dat2[:, 0], dat1[:, 0], dat1[:, 1])
    y = dat2[:, 1]
    rel_err = (yhat - y) / y
    plt.twinx()
    plt.plot(dat2[:, 0], rel_err, "r--")
    plt.ylim(-0.001, 0.001)

plt.show()
