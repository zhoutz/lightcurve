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
  std::unordered_map<int, std::vector<std::pair<int, double>>> spots;

  double get_theta(int i) const { return 0.5 * d_theta + i * d_theta; }

  double get_phi(int i, int j) const {
    if (i % 2 == 0) {
      return j * d_phi;
    } else {
      return 0.5 * d_phi + j * d_phi;
    }
  }

  GridSpots(int n_theta, int n_phi)
      : N_theta(n_theta), N_phi(n_phi), d_theta(pi / N_theta), d_phi(two_pi / N_phi) {}

  double dOmega(int i) const { return d_theta * d_phi * std::sin(get_theta(i)); }

  template <typename F> void init_map(F &&f) {
    for (int i = 0; i < N_theta; ++i) {
      for (int j = 0; j < N_phi; ++j) {
        double t = f(get_theta(i), get_phi(i, j));
        if (t > 0.) {
          spots[i].emplace_back(j, t);
        }
      }
    }
  }

  void print_map() const {
    for (const auto &[i, jt] : spots) {
      std::cout << "theta: " << get_theta(i) << ", dOmega: " << dOmega(i) << "\n";
      for (const auto &[j, t] : jt) {
        std::cout << "  phi: " << get_phi(i, j) << ", value: " << t << "\n";
      }
    }
  }
};
