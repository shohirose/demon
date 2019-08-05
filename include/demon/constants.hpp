#ifndef DEMON_CONSTANTS_HPP
#define DEMON_CONSTANTS_HPP

constexpr int N = 300;         //粒子数
constexpr int n = 9;           //完全平衡木の深さ, = ceil(log2(N))
constexpr auto q = 512 - 300;  // = pow(2,n) - N
constexpr auto p = ((N - q) / 2);

constexpr double a = 0.02;      // radius of disk
constexpr double e = 0.95;      //粒子同士の反発係数
constexpr double e_wall = 1.0;  //壁との反発係数
constexpr double g = 1.0;       //規格化された重力加速度
constexpr double Xmin = -1.0;   //左の壁の位置
constexpr double Xmax = 1.0;    //右の壁の位置
constexpr double Ymin = 0.0;    //底面の位置
constexpr double Ymax = 1.0;    //セルの最高点の位置
constexpr double U = 0.149;     //床面の振動する速度
constexpr double V0 = 1.0;      //初期条件での速度分布の標準偏差
constexpr int N_cell_x = 32;    // x方向のセルの分割数
constexpr int N_cell_y = 12;    // y方向のセルの分割数
constexpr double T = 20.0;      //シミュレーション終了時刻
constexpr double epsilon = 0.000001;

#endif  // DEMON_CONSTANTS_HPP
