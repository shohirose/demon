#include "demon/random.hpp"
#include <algorithm>
#include <functional>

seed_t create_seed() {
  std::random_device rnd;
  seed_t seeds;
  std::generate(seeds.begin(), seeds.end(), std::ref(rnd));
  return seeds;
}

std::mt19937 create_random_device() {
  const auto seeds = create_seed();
  std::seed_seq seq(seeds.begin(), seeds.end());
  return std::mt19937(seq);
}

std::mt19937& random_engine() {
  static thread_local std::mt19937 engine = create_random_device();
  return engine;
}
