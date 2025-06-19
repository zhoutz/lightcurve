#pragma once

#include "linear_interp.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct LensingTable {
  std::vector<double> cos_psi, cos_alpha, lensing_factor, dt;
  LinearInterp cos_psi_interp;

  LensingTable() {
    std::string in_fname = "lensing.txt";
    std::ifstream in_file(in_fname);
    if (!in_file) {
      std::cerr << "Could not open lensing table file: " << in_fname << "\n";
      std::exit(EXIT_FAILURE);
    }
    double cos_psi_, cos_alpha_, lensing_factor_, dt_;
    while (in_file >> cos_psi_ >> cos_alpha_ >> lensing_factor_ >> dt_) {
      cos_psi.push_back(cos_psi_);
      cos_alpha.push_back(cos_alpha_);
      lensing_factor.push_back(lensing_factor_);
      dt.push_back(dt_);
    }
    in_file.close();
    cos_psi.shrink_to_fit();
    cos_alpha.shrink_to_fit();
    lensing_factor.shrink_to_fit();
    dt.shrink_to_fit();

    cos_psi_interp.reset(cos_psi.data(), cos_psi.size());
  }

  auto cal_lens_of_cos_psi(double cos_psi_) {
    int i = cos_psi_interp.hunt(cos_psi_);
    double cos_alpha_ = cos_psi_interp.lin_interp(cos_psi_, cos_psi[i], cos_psi[i + 1],
                                                  cos_alpha[i], cos_alpha[i + 1]);
    double lensing_factor_ = cos_psi_interp.lin_interp(cos_psi_, cos_psi[i], cos_psi[i + 1],
                                                       lensing_factor[i], lensing_factor[i + 1]);
    double dt_ = cos_psi_interp.lin_interp(cos_psi_, cos_psi[i], cos_psi[i + 1], dt[i], dt[i + 1]);
    return std::make_tuple(cos_alpha_, lensing_factor_, dt_);
  }
};
