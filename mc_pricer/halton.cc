/**
 * @file halton.cc
 * @brief Implementation of Halton sequence generator.
 */
#include "halton.h"
#include <cstddef>

namespace pricer {

    /** @brief Internal primality check used to select bases. */
    static bool is_prime(unsigned n){
        if (n < 2) return false;
        if (n % 2 == 0) return n == 2;
        for (unsigned d = 3; d * d <= n; d += 2)
            if (n % d == 0) return false;
        return true;
    }

    /** @brief Construct Halton with prime bases for the first @p dim primes. */
    Halton::Halton(unsigned dim) {
        bases.reserve(dim);
        unsigned candidate = 2;
        while (bases.size() < dim) {
            if (is_prime(candidate)) bases.push_back(candidate);
            ++candidate;
        }
    }

    /** @brief Compute radical inverse in base @p b. */
    double Halton::radical_inverse(unsigned long long n, unsigned b) {
        double inv = 1.0 / b, f = inv, x = 0.0;
        while (n) {
            x += (n % b) * f;
            n /= b;
            f *= inv;
        }
        return x;
    }

    /** @brief Fill @p u with the @p i-th Halton point. */
    void Halton::point(unsigned long long i, std::vector<double>& u) const {
        for (std::size_t d = 0; d < bases.size(); ++d)
            u[d] = radical_inverse(i, bases[d]);
    }

} 