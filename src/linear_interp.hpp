#pragma once

#include <algorithm>

struct LinearInterp {
  int n, jsav;
  double const *xx;

  LinearInterp() : n(0), jsav(0), xx(nullptr) {}

  LinearInterp(double const *x, double const *y, int n_)
      : n(n_), jsav(0), xx(x) {}

  void reset(double const *x, int n_) {
    jsav = n_ - 1;
    n = n_;
    xx = x;
  }

  int hunt(double x) {
    int jl = jsav, jm, ju, inc = 1;
    bool ascnd = (xx[n - 1] >= xx[0]);
    if (jl < 0 || jl > n - 1) {
      jl = 0;
      ju = n - 1;
    } else {
      if (x >= xx[jl] == ascnd) {
        for (;;) {
          ju = jl + inc;
          if (ju >= n - 1) {
            ju = n - 1;
            break;
          } else if (x < xx[ju] == ascnd)
            break;
          else {
            jl = ju;
            inc += inc;
          }
        }
      } else {
        ju = jl;
        for (;;) {
          jl = jl - inc;
          if (jl <= 0) {
            jl = 0;
            break;
          } else if (x >= xx[jl] == ascnd)
            break;
          else {
            ju = jl;
            inc += inc;
          }
        }
      }
    }
    while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if (x >= xx[jm] == ascnd)
        jl = jm;
      else
        ju = jm;
    }
    jsav = jl;
    return std::max(0, std::min(n - 2, jl));
  }

  static double lin_interp(double x, double x1, double x2, double y1,
                           double y2) {
    if (x1 == x2) {
      return y1;
    } else {
      return y1 + ((x - x1) / (x2 - x1)) * (y2 - y1);
    }
  }

  static double log_interp(double x, double x1, double x2, double y1,
                           double y2) {
    if (x1 == x2) {
      return y1;
    } else {
      double k = std::log(y2 / y1) / std::log(x2 / x1);
      return y1 * std::pow(x / x1, k);
    }
  }
};
