#pragma once

#include <vector>

// https://numpy.org/doc/stable/reference/generated/numpy.gradient.html
inline std::vector<double> gradient(std::vector<double> const &y, std::vector<double> const &x) {
  int n = y.size();
  if (n == 1) {
    return std::vector<double>{1};
  }
  std::vector<double> dydx(n);
  dydx[0] = (y[1] - y[0]) / (x[1] - x[0]);
  dydx[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
  for (int i = 1; i <= n - 2; ++i) {
    double hd = x[i + 1] - x[i];
    double hs = x[i] - x[i - 1];
    double hh = hd + hs;
    double wl = -hd / (hs * hh);
    double wm = (hd - hs) / (hs * hd);
    double wr = hs / (hd * hh);
    dydx[i] = wl * y[i - 1] + wm * y[i] + wr * y[i + 1];
  }
  return dydx;
}
