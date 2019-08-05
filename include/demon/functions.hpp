#ifndef DEMON_FUNCTIONS_HPP
#define DEMON_FUNCTIONS_HPP

#include <vector>
#include "demon/cell.hpp"
#include "demon/node.hpp"
#include "demon/particle.hpp"

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

void CBT_build(std::vector<std::vector<Node>>& nodes,
               std::vector<Particle>& particles);

void CBT_update(std::vector<std::vector<Node>>& nodes, double time_new,
                int i_new);

void G1(Particle& particle, int j);
void G2(Particle& p1, Particle& p2);

double NextEvent(std::vector<Particle>& particles,
                 std::vector<std::vector<Cell>>& cells,
                 std::vector<std::vector<Node>>& nodes, int i_current,
                 int j_current);

double t_cell_update(const Particle& particle, int j_current, double t_cell_old,
                     double& v_max);

void MaskUpdate(std::vector<Particle>& particles,
                std::vector<std::vector<Cell>>& cells,
                std::vector<std::vector<Node>>& nodes, int i_current, double t);

double EEPGM(std::vector<Particle>& particles,
             std::vector<std::vector<Cell>>& cells,
             std::vector<std::vector<Node>>& nodes, double t, double& v_max);

#endif  // FUNCTIONS_HPP
