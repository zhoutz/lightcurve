#pragma once

#include "unit.hpp"
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

struct GridSpots {
  int N_theta, N_phi;
  double d_theta, d_phi;
  std::vector<double> theta, phi;
  std::unordered_map<int, std::vector<std::pair<int, double>>> spots;

  GridSpots(int n_theta, int n_phi)
      : N_theta(n_theta), N_phi(n_phi), d_theta(pi / N_theta), d_phi(two_pi / N_phi),
        theta(n_theta), phi(n_phi) {
    for (int i = 0; i < N_theta; ++i) {
      theta[i] = 0.5 * d_theta + i * d_theta;
    }
    for (int j = 0; j < N_phi; ++j) {
      phi[j] = 0.5 * d_phi + j * d_phi;
    }
  }

  double dOmega(int i) const { return d_theta * d_phi * std::sin(theta[i]); }

  template <typename F> void init_map(F &&f) {
    for (int i = 0; i < N_theta; ++i) {
      double domega = dOmega(i);
      for (int j = 0; j < N_phi; ++j) {
        double t = f(theta[i], phi[j]);
        if (t > 0.) {
          spots[i].emplace_back(j, t);
        }
      }
    }
  }

  void print_map() const {
    for (const auto &[i, jt] : spots) {
      std::cout << "theta: " << theta[i] << ", dOmega: " << dOmega(i) << "\n";
      for (const auto &[j, t] : jt) {
        std::cout << "  phi: " << phi[j] << ", value: " << t << "\n";
      }
    }
  }
};
