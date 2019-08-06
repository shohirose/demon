#ifndef DEMON_CONSTANTS_HPP
#define DEMON_CONSTANTS_HPP

constexpr int NP = 300;  // Number of particles
constexpr int BTD = 9;   // Depth of the balanced tree = ceil(log2(NP))
// Parameters for the balanced tree
constexpr auto Q = 512 - 300;  // = pow(2,BTD) - NP
constexpr auto P = (NP - Q) / 2;
constexpr double RDISK = 0.02;  // Radius of disk
constexpr double PREST = 0.95;  // Coefficient of restitution between particles
// Coefficient of restitution between a particle and a wall
constexpr double WREST = 1.0;
constexpr double GRAVITY = 1.0;  // Normalized gravity acceleration
constexpr double XMIN = -1.0;    // X of the left wall
constexpr double XMAX = 1.0;     // X of the right wall
constexpr double YMIN = 0.0;     // Y of the bottom wall
constexpr double YMAX = 1.0;     // Y of the top wall
constexpr double U = 0.149;      // Oscillation vecolity of the bottom wall
// Standard deviation of velocity distribution at the initial state
constexpr double V0 = 1.0;
constexpr int NCELLX = 32;        // Number of cells in x direction
constexpr int NCELLY = 12;        // Number of cells in y direction
constexpr double TEND = 20.0;     // End time of simulation
constexpr double EPS = 0.000001;  // Tolerance

#endif  // DEMON_CONSTANTS_HPP
