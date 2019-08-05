#ifndef DEMON_RANDOM_HPP
#define DEMON_RANDOM_HPP

#include <array>
#include <random>

using seed_t =
    std::array<std::random_device::result_type, std::mt19937::state_size>;

seed_t create_seed();

std::mt19937 create_random_device();

std::mt19937& random_engine();

#endif  // DEMON_RANDOM_HPP
