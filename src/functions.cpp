#include "demon/functions.hpp"
#include <cmath>
#include <iostream>
#include "demon/constants.hpp"
#include "demon/random.hpp"

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
  std::uniform_real_distribution<> udist(0.0, 1.0);
  std::normal_distribution<> ndist(0.0, V0);
  auto& rd = random_engine();

  for (size_t i = 0; i < particles.size(); i++) {
    const auto x = udist(rd);
    particles[i].x = (Xmin + a) * x + (Xmax - a) * (1 - x);
    particles[i].y = (Ymin + a) * x + (0.5 * Ymax - a) * (1 - x);
    while (set(particles, i) == 0) {
      const auto x2 = udist(rd);
      particles[i].x = (Xmin + a) * x2 + (Xmax - a) * (1 - x2);
      const auto x3 = udist(rd);
      particles[i].y = (Ymin + a) * x3 + (Ymax - a) * (1 - x3);
    }
    particles[i].u = ndist(rd);
    particles[i].v = ndist(rd);
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

void CBT_build(std::vector<std::vector<Node>>& nodes,
               std::vector<Particle>& particles) {
  // initialization for bottom nodes
  for (int i = 0; i < 2 * p + 2 * q; i++) {
    if (i < 2 * p) {
      nodes[n][i].time = particles[i].event.time;
      nodes[n][i].number = i;
    } else {
      nodes[n][i].time = particles[2 * p + (i - 2 * p) / 2].event.time;
      nodes[n][i].number = 2 * p + (i - 2 * p) / 2;
    }
  }

  // tournament
  for (int n_index = n - 1; n_index >= 0; n_index--) {
    for (int i = 0; i <= std::pow(2, n_index) - 1; i++) {
      nodes[n_index][i].left = &nodes[n_index + 1][2 * i];
      nodes[n_index + 1][2 * i].parent = &nodes[n_index][i];
      nodes[n_index][i].right = &nodes[n_index + 1][2 * i + 1];
      nodes[n_index + 1][2 * i + 1].parent = &nodes[n_index][i];
      if (nodes[n_index + 1][2 * i].time <=
          nodes[n_index + 1][2 * i + 1].time) {
        nodes[n_index][i].time = nodes[n_index + 1][2 * i].time;
        nodes[n_index][i].number = nodes[n_index + 1][2 * i].number;
      } else {
        nodes[n_index][i].time = nodes[n_index + 1][2 * i + 1].time;
        nodes[n_index][i].number = nodes[n_index + 1][2 * i + 1].number;
      }
    }
  }
}

void CBT_update(std::vector<std::vector<Node>>& entry, double time_new,
                int i_new) {
  Node* entry_now;
  if (i_new < 2 * p) {
    entry[n][i_new].time = time_new;
    entry_now = &entry[n][i_new];
  } else {
    entry[n][2 * i_new - 2 * p].time = time_new;
    entry[n][2 * i_new - 2 * p + 1].time = time_new;
    entry_now = &entry[n][2 * i_new - 2 * p];
  }

  while (entry_now->parent != nullptr) {
    entry_now = entry_now->parent;
    if (entry_now->left->time < entry_now->right->time) {
      entry_now->time = entry_now->left->time;
      entry_now->number = entry_now->left->number;
    } else {
      entry_now->time = entry_now->right->time;
      entry_now->number = entry_now->right->number;
    }
  }
}

void G1(Particle& particle, int j) {
  if ((j == -1) || (j == -2)) {  // collision with R or L wall
    particle.u = -e_wall * particle.u;
    if (j == -1) {
      particle.x = Xmax - a - epsilon;
      //このepsilon処理はgetcell_xのときなどに必要になる
    } else {
      particle.x = Xmin + a + epsilon;
    }
  } else if (j == -3) {  // collision with Bottom wall
    particle.v = (1 + e_wall) * U - e_wall * particle.v;
    particle.y = Ymin + a + epsilon;
  }
}

void G2(Particle& p1, Particle& p2) {
  const auto d = r_distance(p1, p2);
  const auto u1 = p1.u;
  const auto v1 = p1.v;
  const auto u2 = p2.u;
  const auto v2 = p2.v;
  const auto cx = (p1.x - p2.x) / d;
  const auto cy = (p1.y - p2.y) / d;
  p1.u = 0.5 * (1 + e) * ((u2 - u1) * cx + (v2 - v1) * cy) * cx + u1;
  p1.v = 0.5 * (1 + e) * ((u2 - u1) * cx + (v2 - v1) * cy) * cy + v1;
  p2.u = 0.5 * (1 + e) * ((u1 - u2) * cx + (v1 - v2) * cy) * cx + u2;
  p2.v = 0.5 * (1 + e) * ((u1 - u2) * cx + (v1 - v2) * cy) * cy + v2;
}

double NextEvent(std::vector<Particle>& particles,
                 std::vector<std::vector<Cell>>& cells,
                 std::vector<std::vector<Node>>& nodes, int i_current,
                 int j_current) {
  double t = particles[i_current].event.time;
  Free_evolution(particles[i_current],
                 t - particles[i_current].tau);  // i_currentの時間発展
  if (j_current >= 0) {                          // Disk Disk Collision
    Free_evolution(particles[j_current],
                   t - particles[j_current].tau);    // j_currentの時間発展
    G2(particles[i_current], particles[j_current]);  //粒子同士の衝突処理
  }
  if (j_current < 0) {                    // Disk Wall Collision
    G1(particles[i_current], j_current);  //壁との衝突処理
  }
  particles[i_current].event =
      Predictions(particles, cells, t, i_current);  // i_currentのイベント更新
  CBT_update(nodes, particles[i_current].event.time,
             i_current);  // i_currentのnodeアップデート
  if (j_current >= 0) {   // j_currentについても同様
    particles[j_current].event = Predictions(particles, cells, t, j_current);
    CBT_update(nodes, particles[j_current].event.time, j_current);
  }
  return t;
}

double t_cell_update(const Particle& particle, int j_current, double t_cell_old,
                     double& v_max) {
  double t_cell, dt_cell;
  // double cell_length_x = (Xmax - Xmin) / (double)N_cell_x;
  double cell_length_y = (Ymax - Ymin) / (double)N_cell_y;
  if (j_current == -3) {
    const auto v = particle.velocity();
    if (v_max < v) v_max = v;
  }
  dt_cell = (cell_length_y - 2.0 * a) / (2.0 * (v_max));
  t_cell = t_cell_old + dt_cell;
  return t_cell;
}

void MaskUpdate(std::vector<Particle>& particles,
                std::vector<std::vector<Cell>>& cells,
                std::vector<std::vector<Node>>& nodes, int i_current,
                double t) {
  double cell_length_x = (Xmax - Xmin) / (double)N_cell_x,
         cell_length_y = (Ymax - Ymin) / (double)N_cell_y;
  int cell_x = getcell_x(particles[i_current].x, cell_length_x),
      cell_y = getcell_y(particles[i_current].y, cell_length_y), j;
  for (int c1 = -1; c1 <= 1; c1++) {
    for (int c2 = -1; c2 <= 1; c2++) {
      if ((((cell_x + c1 >= 0) && (cell_x + c1 < N_cell_x)) &&
           (cell_y + c2 >= 0)) &&
          (cell_y + c2 < N_cell_y)) {
        j = cells[cell_x + c1][cell_y + c2].first;
        while (j >= 0) {
          if (particles[j].event.number_particle == i_current) {
            particles[j].event = Predictions(particles, cells, t, j);
            CBT_update(nodes, particles[j].event.time, j);
          }
          j = particles[j].next;
        }
      }
    }
  }
}

double EEPGM(std::vector<Particle>& particles,
             std::vector<std::vector<Cell>>& cells,
             std::vector<std::vector<Node>>& nodes, double t, double& v_max) {
  // double cell_length_x = (Xmax - Xmin) / (double)N_cell_x;
  double cell_length_y = (Ymax - Ymin) / (double)N_cell_y;
  double dt_cell, t_cell;

  for (int i = 0; i < N; i++) {  //現在時刻まで時間発展
    Free_evolution(particles[i], t - particles[i].tau);
  }

  cell_register(particles, cells);  //全粒子をセルに登録し直す
  for (int i = 0; i < N; i++) {  //全粒子についてeventを計算し直す
    particles[i].event = Predictions(particles, cells, t, i);
  }
  CBT_build(nodes, particles);  // CBTも最初から構成
  v_max = Vmax(particles);
  dt_cell = (cell_length_y - 2.0 * a) / (2.0 * v_max);
  t_cell = t + dt_cell;
  return t_cell;
}
