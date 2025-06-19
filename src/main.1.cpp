// #include "grid.hpp"
// #include "lensing_table.hpp"
// #include "unit.hpp"
// #include <algorithm>
// #include <cmath>
// #include <iomanip>
// #include <iostream>

// // blackbody spectrum
// // E in keV
// // kT in keV
// double Ibb(double E, double kT) { return h3c2 * (E * E * E) / (std::exp(E / kT) - 1); }

// double cos_psi(double obs_theta, double spot_theta, double spot_phi) {
//   // obs_theta: angle between the observer and the star's axis
//   // spot_theta: angle between the star's axis and the spot
//   // spot_phi: angle between the observer and the spot
//   return std::cos(obs_theta) * std::cos(spot_theta) +
//          std::sin(obs_theta) * std::sin(spot_theta) * std::cos(spot_phi);
// }

// int main() {
//   double mass = 1.4;  // M_sun;
//   double radius = 12; // km
//   double kT = 0.35;   // temperature in keV
//   double spot_angular_radius = 0.01;
//   double spot_theta = 90. * degree;
//   double obs_theta = 90. * degree;
//   double E_obs = 1.0; // keV
//   double D = 0.2 * kpc_in_km;
//   double frequency_nu = 1.0; // Hz
//   int n_time_samples = 400;

//   double u = mass / radius * schwarzschild_radius_of_sun; // compactness
//   // double v = std::sqrt(1. - u);
//   double dS =
//       2. * pi * spot_angular_radius * spot_angular_radius * (1. - std::cos(spot_angular_radius));
//   double D2 = D * D;
  
//   LensingTable lt;

//   for (int i = 0; i < n_time_samples; ++i) {
//     double total_flux = 0.0;

//     double spot_phi = 2. * pi * i / n_time_samples;
//     double c_psi = cos_psi(obs_theta, spot_theta, spot_phi);
//     double y = 1. - c_psi;

//     double x = (1. - u) * y *
//                (1. + u * u * y * y / 112. - e * 0.01 * u * y * (std::log(1 - 0.5 * y) + 0.5 * y));

//     double lensing_factor =
//         (1. - u) *
//         (1. + u * u * y * y * 3. / 112. -
//          0.01 * e * u * y * (2. * std::log(1 - 0.5 * y) + y * (1. - 0.75 * y) / (1 - 0.5 * y)));

//     double cos_alpha = 1. - x;

//     if (cos_alpha < 0) {
//       continue;
//     }

//     total_flux += v * Ibb(E_emit, kT) * cos_alpha * lensing_factor * dS / D2;

//     std::cout << std::setprecision(16);
//     std::cout << total_flux << "\n";
//   }
// }
