#include "lensing_table.hpp"
#include "unit.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

double cal_cos_psi(double obs_theta, double spot_theta, double spot_phi) {
  return std::cos(obs_theta) * std::cos(spot_theta) +
         std::sin(obs_theta) * std::sin(spot_theta) * std::cos(spot_phi);
}

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
  // double u = 0.1723 * 2;
  // std::cout << "u: " << u << "\n";
  double uu = std::sqrt(1. - u);
  double dS = two_pi * Rstar * Rstar * (1. - std::cos(spot_angular_radius));
  double D2 = D * D;
  double frequency_Omega = frequency_nu * two_pi; // rad/s

  LensingTable lt;

  std::string out_fname = "sd1.txt";
  std::ofstream out_file(out_fname);

  std::ofstream dbg_file("dbg.txt");

  for (int i = 0; i < n_phase; ++i) {
    double phase = double(i) / n_phase;
    double spot_phi = two_pi * phase;
    double cos_psi = cal_cos_psi(obs_theta, spot_theta, spot_phi);
    if (cos_psi < lt.cos_psi.back()) continue;
    auto [cos_alpha, lensing_factor, dt] = lt.cal_lens_of_cos_psi(cos_psi);

    double sin_alpha_over_sin_psi =
        cos_psi == 1. ? std::sqrt(lensing_factor)
                      : std::sqrt((1. - cos_alpha * cos_alpha) / (1. - cos_psi * cos_psi));

    double beta = frequency_Omega * Rstar * std::sin(spot_theta) / uu / c_in_km_s;
    double gamma = 1. / std::sqrt(1. - beta * beta);
    double cos_xi = -sin_alpha_over_sin_psi * std::sin(obs_theta) * std::sin(spot_phi);
    double delta = 1. / (gamma * (1. - beta * cos_xi));
    double delta3 = delta * delta * delta;
    double E_emit = E_obs / (delta * uu);
    double delta_phase = dt * frequency_nu;

    double total_flux =
        uu * delta3 * Ibb(E_emit, kT) * cos_alpha * lensing_factor * dS / D2 / gamma;

    dbg_file << std::setprecision(16);
    dbg_file << phase << " ";
    dbg_file << cos_psi << '\n';

    std::cout << std::setprecision(10);
    std::cout << "phase: " << phase << ", ";
    std::cout << "delta_phase: " << delta_phase << ",";
    // std::cout << "cos_psi: " << cos_psi << ", ";
    // std::cout << "cos_alpha: " << cos_alpha << ", ";
    // std::cout << "lensing_factor: " << lensing_factor << ", ";
    // std::cout << "dt: " << dt << ", ";
    // std::cout << "E_emit: " << E_emit << ", ";
    std::cout << "total_flux: " << total_flux << ", ";
    std::cout << "\n";

    out_file << std::setprecision(16);
    out_file << phase << " " << total_flux << " " << delta_phase << "\n";
    // if (i == 2) break; // for testing, remove this line to calculate all phases
  }
}
