/**
 * @file rng.h
 * @brief Simple RNG wrapper providing normal and uniform draws.
 */
#ifndef PRICER_RNG_H
#define PRICER_RNG_H

#include <random>
#include <algorithm>
#include <limits>

namespace pricer {
/**
 * @brief RNG wrapper (mt19937_64) with normal and uniform helpers.
 */
struct RNG {
    std::mt19937_64 gen;
    std::normal_distribution<double> nd{0.0, 1.0};
    std::uniform_real_distribution<double> ud{0.0, 1.0};
    /** @brief Construct RNG with a 64-bit seed. */
    explicit RNG(uint64_t seed) : gen(seed) {}
    /** @brief Draw a standard normal variate. */
    inline double normal() { return nd(gen); }
    /** @brief Draw a uniform variate in [0,1) (capped just below 1). */
    inline double uniform() { return std::min(ud(gen), std::nextafter(1.0, 0.0)); }
};
} 
#endif
