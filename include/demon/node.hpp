#ifndef DEMON_NODE_HPP
#define DEMON_NODE_HPP

struct Node {   //完全平衡木のノードの構造体
  int number;   //対応する粒子番号
  double time;  //対応する粒子の予測された最短衝突時刻
  Node *left;   //上と左右をつなぐ
  Node *right;
  Node *parent;
};

#endif  // DEMON_NODE_HPP
