#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>
#include "cell.hpp"
#include "particle.hpp"

/// @brief Returns uniform distribution
double Uniform();

/// @brief Returns normal distribution
/// @param[in] mu Mean
/// @param[in] sigma Standard deviation
double rand_normal(double mu, double sigma);

int set(std::vector<Particle>& particles, int i);

void status_initialize(std::vector<Particle>& particles);

double Vmax(const std::vector<Particle>& particles);

int getcell_x(double x, double cell_length_x);

int getcell_y(double y, double cell_length_y);

void cell_register(std::vector<Particle>& particles,
                   std::vector<std::vector<Cell>>& cells);

double T_DDC(const Particle& p1, const Particle& p2, double t);

double T_DWC(const Particle& particle, double t, int j);

struct Event Predictions(std::vector<Particle>& particles,
                         std::vector<std::vector<Cell>>& cells, double t,
                         int i);

#endif  // FUNCTIONS_HPP