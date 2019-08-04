#include "particle.hpp"
#include <cmath>
#include "constants.hpp"

double r_distance(const Particle& p1, const Particle& p2) {
  const auto dx = p1.x - p2.x;
  const auto dy = p1.y - p2.y;
  return std::sqrt(dx * dx + dy * dy);
}

double v_distance(const Particle& p1, const Particle& p2) {
  const auto du = p1.u - p2.u;
  const auto dv = p1.v - p2.v;
  return std::sqrt(du * du + dv * dv);
}

void Free_evolution(Particle& p, double t) {
  p.x += (p.u) * t;
  p.y += (p.v) * t - 0.5 * g * t * t;
  p.v += -g * t;
  p.tau += t;  //ここはうまくいっているか確認が必要
}