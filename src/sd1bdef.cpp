#include "gradient.hpp"
#include "grid.hpp"
#include "lensing_table.hpp"
#include "linear_interp.hpp"
#include "unit.hpp"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

double frequency_nu, spot_theta, obs_theta;

double temperature_profile(double grid_theta, double grid_phi) {
  double kT = 0.35; // temperature in keV
  double spot_angular_radius = 1.0;
  // double spot_theta = 20. * degree;
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

void fill_phase_o(std::vector<double> &phase_o) {
  int n = phase_o.size();
  double phase_min = 0, phase_max = double(n - 1) / n, dphase = 1. / n;

  auto do_fill = [&](int beg, int end) {
    double phase0, phase1;

    if (beg == 0)
      phase0 = phase_min;
    else
      phase0 = phase_o[beg - 1] + dphase;

    if (end == n - 1)
      phase1 = phase_max;
    else
      phase1 = phase_o[end + 1] - dphase;

    for (int i = beg; i <= end; ++i) {
      phase_o[i] = phase0 + (phase1 - phase0) * (double(i - beg) / (end - beg));
    }
  };

  int beg = -1;
  for (int i = 0; i < n; ++i) {
    if (phase_o[i] == -1) {
      if (beg == -1) beg = i;
    } else {
      if (beg != -1) {
        int end = i - 1;
        do_fill(beg, end);
        beg = -1;
      }
    }
  }
  if (beg != -1) {
    int end = n - 1;
    do_fill(beg, end);
  }
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " <frequency_nu> <spot_theta> <obs_theta>\n";
  } else {
    frequency_nu = std::atof(argv[1]);
    spot_theta = std::atof(argv[2]) * degree;
    obs_theta = std::atof(argv[3]) * degree;
  }

  double Mstar = 1.4; // M_sun;
  double Rstar = 12;  // km
  double E_obs = 1.0; // keV

  // double obs_theta = 80. * degree;
  double D = 0.2 * kpc_in_km;
  // double frequency_nu = 400.0; // Hz
  int n_phase = 250;

  double u = Mstar / Rstar * schwarzschild_radius_of_sun; // compactness
  double uu = std::sqrt(1. - u);
  double R2 = Rstar * Rstar;
  double D2 = D * D;
  double frequency_Omega = frequency_nu * two_pi; // rad/s

  LensingTable lt;
  GridSpots grid_spots(250, n_phase);
  grid_spots.init_map(temperature_profile);

  std::vector<double> total_fluxes(n_phase, 0.0), phase_output(n_phase);
  std::vector<double> fluxes_div_I(n_phase), redshift_factors(n_phase);
  std::vector<double> phase_s(n_phase), phase_o(n_phase);
  std::vector<double> single_patch_fluxes(n_phase);
  LinearInterp li;
  li.reset(phase_o.data(), phase_o.size());

  for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
    phase_output[i_phase] = double(i_phase) / n_phase;
  }

  auto time_start = std::chrono::high_resolution_clock::now();

  for (auto const &[i, jt] : grid_spots.spots) {
    double spot_theta = grid_spots.get_theta(i);
    double dOmega = grid_spots.dOmega(i);
    double dS = dOmega * R2;

    double cc = std::cos(obs_theta) * std::cos(spot_theta);
    double ss = std::sin(obs_theta) * std::sin(spot_theta);

    for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
      double phase = double(i_phase) / n_phase;
      double spot_phi = two_pi * phase + grid_spots.get_phi(i, 0);
      double cos_psi = cc + ss * std::cos(spot_phi);
      if (cos_psi < lt.cos_psi.back()) {
        fluxes_div_I[i_phase] = 0;
        redshift_factors[i_phase] = 0;
        phase_s[i_phase] = phase;
        phase_o[i_phase] = -1;
      } else {
        auto [cos_alpha, lensing_factor, dt] = lt.cal_lens_of_cos_psi(cos_psi);
        double sin_alpha_over_sin_psi =
            cos_psi == 1. ? std::sqrt(lensing_factor)
                          : std::sqrt((1. - cos_alpha * cos_alpha) / (1. - cos_psi * cos_psi));
        double beta = frequency_Omega * Rstar * std::sin(spot_theta) / uu / c_in_km_s;
        double gamma = 1. / std::sqrt(1. - beta * beta);
        double cos_xi = -sin_alpha_over_sin_psi * std::sin(obs_theta) * std::sin(spot_phi);
        double delta = 1. / (gamma * (1. - beta * cos_xi));
        double delta3 = delta * delta * delta;
        double delta_phase = dt * frequency_nu;

        fluxes_div_I[i_phase] = uu * delta3 * cos_alpha * lensing_factor * dS / D2 / gamma;
        redshift_factors[i_phase] = 1. / (delta * uu);
        phase_s[i_phase] = phase;
        phase_o[i_phase] = phase + delta_phase;
      }
    }

    fill_phase_o(phase_o);
    auto dphase = gradient(phase_s, phase_o);

    for (auto const &[j, kt] : jt) {
      for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
        if (fluxes_div_I[i_phase] == 0) {
          single_patch_fluxes[i_phase] = 0;
        } else {
          double target_flux = Ibb(E_obs * redshift_factors[i_phase], kt) * fluxes_div_I[i_phase];
          single_patch_fluxes[i_phase] = target_flux * dphase[i_phase];
        }
      }

      for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
        double phase_shift = double(j) / n_phase;
        double patch_phase = std::fmod(phase_output[i_phase] + phase_shift + 1., 1.);
        int i = li.hunt(patch_phase);
        double flux_output = li.lin_interp(patch_phase, phase_o[i], phase_o[i + 1],
                                           single_patch_fluxes[i], single_patch_fluxes[i + 1]);
        total_fluxes[i_phase] += flux_output;
      }
    }
  }

  auto time_end = std::chrono::high_resolution_clock::now();
  double t_ms = std::chrono::duration<double, std::milli>(time_end - time_start).count();
  std::cout << "Time taken: " << t_ms << " ms\n";

  std::ofstream out_file("sd1.txt");
  out_file << std::setprecision(16);
  for (int i_phase = 0; i_phase < n_phase; ++i_phase) {
    out_file << phase_output[i_phase] << " " << total_fluxes[i_phase] << "\n";
  }
}
