#include <iostream>
#include "constants.hpp"
#include "functions.hpp"
#include "node.hpp"

int main(void) {
  FILE *fp_position, *fp_height;  //ファイルの生成
  char name_position[256];
  sprintf(name_position, "position(N=%d).txt",
          N);  //一定時刻ごとに粒子の位置を保存
  if ((fp_position = fopen(name_position, "w")) == NULL) {
    printf("file_check open error\n");
  }
  if ((fp_height = fopen("height.txt", "w")) ==
      NULL) {  //粒子の高さの平均値を記録
    printf("file open error\n");
  }
  int i_current;
  int j_current;  //現在注目している粒子のペア,j_current<0:壁,j_current>=0:粒子
  double v_max = 0.0;  //最大速度を保存、セルの更新のために必要
  // double cell_length_x = (Xmax - Xmin) / (double)N_cell_x;
  double cell_length_y = (Ymax - Ymin) / (double)N_cell_y;
  double t = 0.0;
  // double dt = 0.01;
  double trec = 0.0;
  double dtrec = (double)T / 200.0;  // 200枚の画像からgifを生成
  double t_cell = 0.;
  double t_cell_old = 0.0;  //セルの更新時刻
  double height;

  std::vector<Particle> particles(N);

  std::vector<std::vector<Cell>> cells(N_cell_x);
  for (auto &&ci : cells) ci.resize(N_cell_y);

  status_initialize(particles);  //位置や速度の初期化

  v_max = Vmax(particles);
  t_cell = (cell_length_y - 2.0 * a) /
           (2.0 * v_max);  //この時間までにマスク外からの衝突はありえない
  cell_register(particles, cells);  //粒子をセルに登録する,nextofの初期化

  std::vector<std::vector<Node>> nodes(n + 1);
  for (auto &&ni : nodes) ni.resize(2 * p + 2 * q);
  //完全平衡木(あるいはトーナメント)を表す構造体

  for (int i = 0; i < N; i++) {
    particles[i].event = Predictions(particles, cells, t,
                                     i);  //それぞれの粒子の次のイベントを予測
  }

  CBT_build(nodes, particles);  // Complete Binary Treeを組み立てる
  printf("set up ok\n");

  while (t <= T) {
    // NEXT EVENTの検索
    i_current =
        nodes[0][0].number;  //決勝のノードは最短の時間で衝突する粒子を示す
    j_current =
        particles[i_current]
            .event
            .number_particle;  // i_currentの衝突相手(これは壁の可能性もある)
    t = NextEvent(particles, cells, nodes, i_current,
                  j_current);  // NEXT EVENTを処理しtとparticle,cell,nodeを更新
    t_cell = t_cell_update(particles[i_current], j_current, t_cell_old,
                           v_max);  // t_cellとv_maxの更新

    // i_current,j_currentと衝突する粒子がいた場合はその粒子のeventはinvalidになってしまうので新しくeventを作る
    //そのような粒子は同じマスク内にしか存在しないはずなのでその中で探索
    MaskUpdate(particles, cells, nodes, i_current,
               t);  // i_currentの周りの粒子でinvalidなものがあればアップデート
    if (j_current >= 0) {  // jについても同様
      MaskUpdate(particles, cells, nodes, j_current, t);
    }

    // EEPGM マスク外の粒子とも衝突する可能性が生じるので登録し直す
    if (t >= t_cell) {
      t_cell_old = t;
      t_cell = EEPGM(particles, cells, nodes, t, v_max);
      //床に粒子がめり込んでいたらこのエラーが生じる
      for (int i = 0; i < N; i++) {
        if (particles[i].y < Ymin + a - epsilon) {
          printf("i=%d:error\n", i);
          printf("%lf %lf %lf %lf\n", particles[i].x, particles[i].y,
                 particles[i].u, particles[i].v);
          printf("%lf %d %d\n", particles[i].event.time,
                 particles[i].event.number_particle,
                 particles[i].event.number_col);
          G1(particles[i], -3);
          particles[i].event = Predictions(particles, cells, t, i);
          CBT_update(nodes, particles[i].event.time, i);
          MaskUpdate(particles, cells, nodes, i, t);
        }
      }
    }
    //粒子の位置の出力
    if ((t > trec) && (t < T)) {
      t_cell_old = t;
      t_cell = EEPGM(particles, cells, nodes, t, v_max);
      printf("t = %lf, v_max = %lf\n", t, v_max);
      height = 0.0;
      for (int i = 0; i < N; i++) {
        fprintf(fp_position, "%lf %lf\n", particles[i].x, particles[i].y);
        height += particles[i].y / (double)N;
      }
      fprintf(fp_position, "\n\n");
      fprintf(fp_height, "%lf %lf\n", t, height);
      trec += dtrec;
    }
  }

  fclose(fp_position);
  fclose(fp_height);
  return 0;
}