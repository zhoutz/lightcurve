#pragma once

#include <vector>
#include <fstream>

struct Spot{
  double theta;
  double phi;
};

inline
std::vector<Spot> read_grid(const std::string& filename) {
  std::vector<Spot> grid;
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  double theta, phi;
  while (file >> theta >> phi) {
    grid.push_back({theta, phi});
  }

  file.close();
  return grid;
}
