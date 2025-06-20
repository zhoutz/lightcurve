import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1, 1, 100)
y = np.arctan(x)

dydx = np.gradient(y, x)
dydx_ = 1 / (x**2 + 1)

plt.plot(x, dydx)
plt.plot(x, dydx_)

plt.show()
