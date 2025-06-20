#include "gradient.hpp"
#include "lensing_table.hpp"
#include "unit.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
  double Mstar = 1.4; // M_sun;
  double Rstar = 12;  // km
  double kT = 0.35;   // temperature in keV
  double E_obs = 1.0; // keV
  double spot_angular_radius = 0.01;
  double spot_theta = 90. * degree;
  double obs_theta = 90. * degree;
  double D = 0.2 * kpc_in_km;
  double frequency_nu = 200.0; // Hz
  int n_phase = 800;

  double u = Mstar / Rstar * schwarzschild_radius_of_sun; // compactness
  double uu = std::sqrt(1. - u);
  double dS = two_pi * Rstar * Rstar * (1. - std::cos(spot_angular_radius));
  double D2 = D * D;
  double frequency_Omega = frequency_nu * two_pi; // rad/s
  double cc = std::cos(obs_theta) * std::cos(spot_theta);
  double ss = std::sin(obs_theta) * std::sin(spot_theta);

  LensingTable lt;

  std::vector<double> phase_s(n_phase), phase_o(n_phase), total_fluxes(n_phase);

  for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
    double phase = double(i_phase) / n_phase;
    double spot_phi = two_pi * phase;
    double cos_psi = cc + ss * std::cos(spot_phi);
    if (cos_psi < lt.cos_psi.back()) continue;
    auto [cos_alpha, lensing_factor, dt] = lt.cal_lens_of_cos_psi(cos_psi);

    double sin_alpha_over_sin_psi =
        cos_psi == 1. ? std::sqrt(lensing_factor)
                      : std::sqrt((1. - cos_alpha * cos_alpha) / (1. - cos_psi * cos_psi));

    double beta = frequency_Omega * Rstar * std::sin(spot_theta) / uu / c_in_km_s;
    double gamma = 1. / std::sqrt(1. - beta * beta);
    double cos_xi = -sin_alpha_over_sin_psi * std::sin(obs_theta) * std::sin(spot_phi);
    double delta = 1. / (gamma * (1. - beta * cos_xi));
    double delta2 = delta * delta;
    double delta3 = delta2 * delta;
    double delta4 = delta2 * delta2;
    double E_emit = E_obs / (delta * uu);
    double delta_phase = dt * frequency_nu;

    double total_flux =
        uu * delta3 * Ibb(E_emit, kT) * cos_alpha * lensing_factor * dS / D2 / gamma;

    phase_s[i_phase] = phase;
    phase_o[i_phase] = phase + delta_phase;
    total_fluxes[i_phase] = total_flux;

    std::cout << std::setprecision(10);
    std::cout << "phase: " << phase << ", ";
    std::cout << "delta_phase: " << delta_phase << ",";
    std::cout << "total_flux: " << total_flux << ", ";
    std::cout << "\n";
  }

  auto dphase = gradient(phase_s, phase_o);

  std::ofstream out_file("sd1.txt");
  out_file << std::setprecision(16);
  for (int i = 0; i < n_phase; ++i) {
    out_file << phase_o[i] << " ";
    out_file << total_fluxes[i] * dphase[i] << "\n";
  }
}
