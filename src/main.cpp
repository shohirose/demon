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
  int i_current,
      j_current;  //現在注目している粒子のペア,j_current<0:壁,j_current>=0:粒子
  double v_max = 0.0;  //最大速度を保存、セルの更新のために必要
  double cell_length_x = (Xmax - Xmin) / (double)N_cell_x,
         cell_length_y = (Ymax - Ymin) / (double)N_cell_y;
  double t = 0.0, dt = 0.01, trec = 0.0,
         dtrec = (double)T / 200.0;       // 200枚の画像からgifを生成
  double t_cell = 0.0, t_cell_old = 0.0;  //セルの更新時刻
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
}