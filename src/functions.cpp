#include "functions.hpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include "constants.hpp"

double Uniform() { return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0); }

double rand_normal(double mu, double sigma) {
  double z =
      std::sqrt(-2.0 * std::log(Uniform())) * std::sin(2.0 * M_PI * Uniform());
  return mu + sigma * z;
}

int set(std::vector<Particle>& particles,
        int i) {  // setに成功していれば1,失敗していれば0を返す
  int r = 1;
  double d;

  if (std::fabs(particles[i].x) < a) {
    r = 0;
  }
  for (int j = 1; j <= i - 1; j++) {
    d = r_distance(particles[i], particles[j]);
    if (d <= 2.0 * a) {
      r = 0;
      break;
    }
  }
  return r;
}

void status_initialize(std::vector<Particle>& particles) {
  for (size_t i = 0; i < particles.size(); i++) {
    const auto prob = Uniform();
    particles[i].x = (Xmin + a) * prob + (Xmax - a) * (1 - prob);
    particles[i].y = (Ymin + a) * prob + (0.5 * Ymax - a) * (1 - prob);
    while (set(particles, i) == 0) {
      const auto prob2 = Uniform();
      particles[i].x = (Xmin + a) * prob2 + (Xmax - a) * (1 - prob2);
      const auto prob3 = Uniform();
      particles[i].y = (Ymin + a) * prob3 + (Ymax - a) * (1 - prob3);
    }
    particles[i].u = rand_normal(0.0, V0);
    particles[i].v = rand_normal(0.0, V0);
    particles[i].next = -1;
    particles[i].tau = 0.0;
    particles[i].event.number_col = 0;
    particles[i].event.time = 2.0 * T;
    particles[i].event.number_particle = -1;
  }
}

double Vmax(const std::vector<Particle>& particles) {
  double v_max = 0.0;
  for (auto&& p : particles) {
    const auto v = p.velocity();
    if (v_max < v) v_max = v;
  }
  return v_max;
}

int getcell_x(double x, double cell_length_x) {
  if ((x < Xmin + a) || (Xmax - a < x)) {
    printf("x is out of range\n");
  }
  return (int)((x - Xmin) / cell_length_x);
}

int getcell_y(double y, double cell_length_y) {
  if (y < Ymin) {
    printf("error:y<0(%lf)\n", y);
    return 0;
  } else if (y > Ymax) {
    return N_cell_y - 1;  // Ymaxよりも高い位置の粒子は一番高いセルに登録
  } else {
    return (int)(y / cell_length_y);
  }
}

void cell_register(std::vector<Particle>& particles,
                   std::vector<std::vector<Cell>>& cells) {
  int i, cell_x, cell_y, lastPrev;
  // initialize particle.next and cell
  for (auto&& p : particles) p.next = -1;

  for (cell_x = 0; cell_x < N_cell_x; cell_x++) {
    for (cell_y = 0; cell_y < N_cell_y; cell_y++) {
      cells[cell_x][cell_y].first = -1;
      cells[cell_x][cell_y].last = -1;
    }
  }

  double cell_length_x = (Xmax - Xmin) / (double)N_cell_x,
         cell_length_y = (Ymax - Ymin) / (double)N_cell_y;
  for (i = 0; i < N; i++) {
    cell_x = getcell_x(particles[i].x, cell_length_x);
    cell_y = getcell_y(particles[i].y, cell_length_y);
    lastPrev = cells[cell_x][cell_y].last;
    cells[cell_x][cell_y].last = i;

    if (lastPrev == -1) {
      cells[cell_x][cell_y].first = i;
    } else {
      particles[lastPrev].next = i;
    }
  }
}

double T_DDC(const Particle& p1, const Particle& p2, double t) {
  double r_relative, v_relative, b, hoge;
  double tau = t;
  double x1 = p1.x, x2 = p2.x, y1 = p1.y, y2 = p2.y;
  double u1 = p1.u, u2 = p2.u, v1 = p1.v, v2 = p2.v;
  r_relative = r_distance(p1, p2);
  v_relative = v_distance(p1, p2);
  b = (x1 - x2) * (u1 - u2) + (y1 - y2) * (v1 - v2);
  hoge =
      b * b - v_relative * v_relative * (r_relative * r_relative - 4.0 * a * a);
  if (hoge > 0) {
    tau += -(b + sqrt(hoge)) / (v_relative * v_relative);
  } else {
    tau += T;
  }
  return tau;
}

double T_DWC(const Particle& particle, double t, int j) {
  double tau = t;
  if (j == -1) {  // collision with RIGHT wall(-1)
    if (particle.u > 0.0) {
      tau += (Xmax - a - particle.x) / particle.u;
    } else {
      tau += 2.0 * T;
    }
  } else if (j == -2) {  // collision with LEFT wall(-2)
    if (particle.u < 0.0) {
      tau += (Xmin + a - particle.x) / particle.u;
    } else {
      tau += 2.0 * T;
    }
  } else if (j == -3) {  // collision with BOTTOM wall(-3)
    tau += (particle.v +
            sqrt(particle.v * particle.v + 2 * g * (particle.y - a))) /
           g;
  }
  if (tau < t) {
    return 2.0 * T;
  } else {
    return tau;
  }
}

struct Event Predictions(std::vector<Particle>& particles,
                         std::vector<std::vector<Cell>>& cells, double t,
                         int i) {
  double t_min = 2.0 * T, t_temp;
  double cell_length_x = (Xmax - Xmin) / (double)N_cell_x,
         cell_length_y = (Ymax - Ymin) / (double)N_cell_y;
  int j_col, j;
  struct Event L;

  for (j = -3; j < 0; j++) {
    t_temp = T_DWC(particles[i], particles[i].tau, j);
    if ((t_temp > t) & (t_temp < t_min)) {
      t_min = t_temp;
      j_col = j;
    }
  }

  int cell_x = getcell_x(particles[i].x, cell_length_x),
      cell_y = getcell_y(particles[i].y, cell_length_y);
  for (int c1 = -1; c1 < 1; c1++) {
    for (int c2 = -1; c2 <= 1; c2++) {
      if ((((cell_x + c1 >= 0) && (cell_x + c1 <= N_cell_x - 1)) &&
           (cell_y + c2 >= 0)) &&
          (cell_y + c2 <= N_cell_y - 1)) {
        j = cells[cell_x + c1][cell_y + c2].first;
        while (j >= 0) {
          Free_evolution(particles[j], particles[i].tau - particles[j].tau);
          t_temp = T_DDC(particles[i], particles[j], particles[i].tau);
          if ((t_temp > t) && (t_temp < t_min)) {
            t_min = t_temp;
            j_col = j;
          }
          j = particles[j].next;
        }
      }
    }
  }
  L.time = t_min;
  L.number_particle = j_col;
  L.number_col = 0;  //今後修正が必要になるかもしれない
  return L;
}