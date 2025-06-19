import numpy as np
from scipy.integrate import quad_vec
import matplotlib.pyplot as plt


def f(x, u, alpha):
    s2 = np.sin(alpha) ** 2
    c2 = 1 - s2
    x2 = x**2
    q = (2 - x2 - u * (1 - x2) ** 2 / (1 - u)) * s2
    tmp1 = np.sqrt(c2 + x2 * q)
    ret1 = x / tmp1
    ret2 = x / (tmp1 * (1 + tmp1))
    return np.array([ret1, ret2])


@np.vectorize
def cal_psi_dt(u, alpha):
    res, err = quad_vec(f, 0, 1, args=(u, alpha))
    s = np.sin(alpha)
    ret1 = 2 * s / np.sqrt(1 - u) * res[0]
    ret2 = 2 * s**2 / (1 - u) * res[1]
    return ret1, ret2


M = 1.4
R = 12
GMsun_c2km = 1.476625038507424
# u = 0.1723*2
u = 2 * GMsun_c2km * M / R
print(f"u = {u:.4f} (M = {M}, R = {R})")
c = 299792458

alpha = np.linspace(0, np.pi / 2, 1024)
psi_dt = cal_psi_dt(u, alpha)
psi = psi_dt[0]
dt = psi_dt[1] * R * 1e3 / c

cos_alpha = np.cos(alpha)
cos_psi = np.cos(psi)
lensing_factor = np.gradient(cos_alpha, cos_psi)

if True:
    out_fname = "lensing.txt"
    np.savetxt(out_fname, np.vstack((cos_psi, cos_alpha, lensing_factor, dt)).T)


data1 = np.loadtxt("apjlab5968/deflection_angle_lensing_factor_test_CU.txt").T
data2 = np.loadtxt("apjlab5968/travel_time_delay_test_CU.txt").T

plt.figure(figsize=(15, 5))
plt.subplot(1, 3, 1)
plt.plot(alpha, psi)
plt.plot(data1[1], data1[0], ".", markersize=1)
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\psi$")

plt.subplot(1, 3, 2)
plt.plot(alpha, lensing_factor)
plt.plot(data1[1], data1[2], ".", markersize=1)
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$d \cos \alpha / d \cos \psi$")

plt.subplot(1, 3, 3)
plt.plot(alpha, dt)
plt.plot(data2[1], data2[2], ".", markersize=1)
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\Delta t$ (s)")

plt.tight_layout()
plt.show()


# err = np.abs(cal_psi(u, data[1]) - data[0])
# rel_err = err / data[0]
# print(f"Max absolute error: {np.max(err):.4e}")


# err = np.abs(psi(u, data[1]) - data[0])
# rel_err = err / data[0]
# print(f"Max absolute error: {np.max(err):.4e}")


plt.show()
