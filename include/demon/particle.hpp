#ifndef DEMON_PARTICLE_HPP
#define DEMON_PARTICLE_HPP

#include <cmath>
#include "demon/event.hpp"

struct Particle {  //粒子に関する情報をまとめる
  double x;        // x position
  double y;        // y position
  double u;        // x velocity
  double v;        // y velocity
  double tau;      //粒子の固有時間を記録，Delayed State
                   // Algorithm(DSA)による高速化のために必要
  double next;     //自分の次に同じセルに入った粒子の番号
  Event event;     //次に予定されるイベント

  auto velocity() const noexcept { return std::sqrt(u * u + v * v); }
};

/// @brief Returns the distance between two particles
/// @param[in] p1 Particle 1
/// @param[in] p2 Particle 2
double r_distance(const Particle& p1, const Particle& p2);

/// @brief Returns the Relative velocity of two particles
/// @param[in] p1 Particle 1
/// @param[in] p2 Particle 2
double v_distance(const Particle& p1, const Particle& p2);

void Free_evolution(Particle& particle, double t);

#endif  // DEMON_PARTICLE_HPP
