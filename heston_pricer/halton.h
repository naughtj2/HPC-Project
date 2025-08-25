/**
 * @file halton.h
 * @brief Halton low-discrepancy sequence generator.
 */
#ifndef PRICER_HALTON_H
#define PRICER_HALTON_H

#include <vector>

namespace pricer {
    /**
     * @brief Generator for Halton sequences using primes as bases.
     */
    struct Halton {
        std::vector<unsigned> bases;
        /** @brief Construct for @p dim dimensions. */
        explicit Halton(unsigned dim);
        /** @brief Radical inverse of index @p n in base @p b. */
        static double radical_inverse(unsigned long long n, unsigned b);
        /** @brief Write the @p i-th point into @p u. */
        void point(unsigned long long i, std::vector<double>& u) const;
    };
} 
#endif
