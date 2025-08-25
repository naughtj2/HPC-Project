/**
 * @file halton.h
 * @brief Halton low-discrepancy sequence generator (base selection and sampling).
 */
#ifndef PRICER_HALTON_H
#define PRICER_HALTON_H

#include <vector>

namespace pricer {

    /**
 * @brief Generator for Halton sequences with prime bases per dimension.
 */
    struct Halton {
        std::vector<unsigned> bases;
        /** @brief Construct with @p dim dimensions (primes chosen as bases). */
        explicit Halton(unsigned dim);
        /** @brief Compute radical inverse of index @p n in base @p b. */
        static double radical_inverse(unsigned long long n, unsigned b);
        /** @brief Generate the @p i-th point into @p u (size must match dimension). */
        void point(unsigned long long i, std::vector<double>& u) const;
    };

} 

#endif
