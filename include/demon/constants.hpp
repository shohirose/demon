#ifndef DEMON_CONSTANTS_HPP
#define DEMON_CONSTANTS_HPP

constexpr int NP = 300;        //粒子数
constexpr int BTD = 9;         //完全平衡木の深さ, = ceil(log2(NP))
constexpr auto Q = 512 - 300;  // = pow(2,D) - NP
constexpr auto P = ((NP - Q) / 2);

constexpr double RDISK = 0.02;   // radius of disk
constexpr double PREST = 0.95;   //粒子同士の反発係数
constexpr double WREST = 1.0;    //壁との反発係数
constexpr double GRAVITY = 1.0;  //規格化された重力加速度
constexpr double XMIN = -1.0;    //左の壁の位置
constexpr double XMAX = 1.0;     //右の壁の位置
constexpr double YMIN = 0.0;     //底面の位置
constexpr double YMAX = 1.0;     //セルの最高点の位置
constexpr double U = 0.149;      //床面の振動する速度
constexpr double V0 = 1.0;       //初期条件での速度分布の標準偏差
constexpr int NCELLX = 32;       // x方向のセルの分割数
constexpr int NCELLY = 12;       // y方向のセルの分割数
constexpr double TEND = 20.0;    //シミュレーション終了時刻
constexpr double EPS = 0.000001;

#endif  // DEMON_CONSTANTS_HPP
