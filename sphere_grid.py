import numpy as np


N_theta = 50
N_phi = 50

theta_min = 0
theta_max = np.pi

phi_min = 0
phi_max = 2 * np.pi

d_theta = (theta_max - theta_min) / N_theta
d_phi = (phi_max - phi_min) / N_phi

theta_grid = np.linspace(theta_min + d_theta / 2, theta_max - d_theta / 2, N_theta)
phi_grid = np.linspace(phi_min + d_phi / 2, phi_max - d_phi / 2, N_phi)

d_area = np.sin(theta_grid) @ np.ones((N_theta, N_phi)) * d_theta * d_phi
total_area = np.sum(d_area)
rel_error = np.abs(total_area - 4 * np.pi) / (4 * np.pi)
print(f"Relative error: {rel_error:.4e}")
