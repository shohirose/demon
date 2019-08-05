#ifndef DEMON_EVENT_HPP
#define DEMON_EVENT_HPP

struct Event {  //ある粒子のイベントの詳細(衝突時刻・相手)を記録
  double time;
  int number_particle;
  int number_col;  //今回は使わない
};

#endif  // DEMON_EVENT_HPP
