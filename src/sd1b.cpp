#include "grid.hpp"
#include "lensing_table.hpp"
#include "unit.hpp"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

double cal_cos_psi(double obs_theta, double spot_theta, double spot_phi) {
  // obs_theta: angle between the observer and the star's axis
  // spot_theta: angle between the star's axis and the spot
  // spot_phi: angle between the observer and the spot
  return std::cos(obs_theta) * std::cos(spot_theta) +
         std::sin(obs_theta) * std::sin(spot_theta) * std::cos(spot_phi);
}

double temperature_profile(double grid_theta, double grid_phi) {
  double kT = 0.35; // temperature in keV
  double spot_angular_radius = 1.0;
  double spot_theta = 90. * degree;
  double spot_phi = 0;

  double s1 = std::sin(grid_theta);
  double c1 = std::cos(grid_theta);
  double s2 = std::sin(spot_theta);
  double c2 = std::cos(spot_theta);

  double cos_rho = s1 * s2 * std::cos(grid_phi - spot_phi) + c1 * c2;
  double cos_angular_radius = std::cos(spot_angular_radius);
  if (cos_rho > cos_angular_radius) {
    return kT;
  } else {
    return -1;
  }
}

int main() {
  double Mstar = 1.4; // M_sun;
  double Rstar = 12;  // km
  double E_obs = 1.0; // keV

  double obs_theta = 90. * degree;
  double D = 0.2 * kpc_in_km;
  double frequency_nu = 1.0; // Hz
  int n_phase = 400;

  double u = Mstar / Rstar * schwarzschild_radius_of_sun; // compactness
  double uu = std::sqrt(1. - u);
  double R2 = Rstar * Rstar;
  double D2 = D * D;
  double frequency_Omega = frequency_nu * two_pi; // rad/s

  LensingTable lt;
  GridSpots grid_spots(100, n_phase);
  grid_spots.init_map(temperature_profile);
  // grid_spots.print_map();
  // std::exit(0);

  std::vector<double> total_fluxes(n_phase, 0.0);

  std::vector<double> fluxes_div_I, redshift_factors;
  for (auto const &[i, jt] : grid_spots.spots) {
    double spot_theta = grid_spots.theta[i];
    double dOmega = grid_spots.dOmega(i);
    double dS = dOmega * R2;

    fluxes_div_I.resize(n_phase, 0.0);
    redshift_factors.resize(n_phase, 0.0);
    for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
      double phase = double(i_phase) / n_phase;
      double spot_phi = two_pi * phase + grid_spots.phi[0];

      double cos_psi = cal_cos_psi(obs_theta, spot_theta, spot_phi);
      if (cos_psi < lt.cos_psi.back()) continue;
      auto [cos_alpha, lensing_factor, dt] = lt.cal_lens_of_cos_psi(cos_psi);

      double sin_alpha_over_sin_psi =
          cos_psi == 1. ? std::sqrt(lensing_factor)
                        : std::sqrt((1. - cos_alpha * cos_alpha) / (1. - cos_psi * cos_psi));

      double beta = frequency_Omega * Rstar * std::sin(spot_theta) / c_in_km_s;
      double gamma = 1. / std::sqrt(1. - beta * beta);
      double cos_xi = -sin_alpha_over_sin_psi * std::sin(obs_theta) * std::sin(spot_phi);
      double delta = 1. / (gamma * (1. - beta * cos_xi));
      double delta2 = delta * delta;
      double delta4 = delta2 * delta2;
      // double E_emit = E_obs / (delta * uu);

      double delta_phase = dt * frequency_nu;
      double delta4gamma = delta4 * gamma;

      fluxes_div_I[i_phase] = uu * delta4gamma * cos_alpha * lensing_factor * dS / D2;
      redshift_factors[i_phase] = 1. / (delta * uu);
    }

    for (auto const &[j, kt] : jt) {
      for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
        if (fluxes_div_I[i_phase] == 0) continue; // Skip if fluxes_div_I is zero
        int target_phase = (i_phase + j) % n_phase;
        double target_flux = Ibb(E_obs * redshift_factors[i_phase], kt) * fluxes_div_I[i_phase];
        total_fluxes[target_phase] += target_flux;
      }
    }
  }

  std::ofstream out_file("sd1.txt");
  out_file << std::setprecision(16);
  for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
    double phase = double(i_phase) / n_phase;
    out_file << phase << " " << total_fluxes[i_phase] << "\n";
  }
}
