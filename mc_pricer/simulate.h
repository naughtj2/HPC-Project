/**
 * @file simulate.h
 * @brief Path simulation under GBM with plain and AAD-aware variants.
 */
#ifndef PRICER_SIMULATE_H
#define PRICER_SIMULATE_H

#include <vector>
#include <type_traits>
#include "options.h"
#include "aad.h"

namespace pricer {

    /** @brief Simulate a GBM path; if T==ADouble and @p tape provided, parameters are active variables.
     *  @tparam T double or ADouble
     *  @param opt Option spec (T, sigma, r, q, steps)
     *  @param Z   Standard normal shocks (size = steps)
     *  @param tape Optional tape to attach active parameters to
     */
    template <typename T>
    std::vector<T> simulate_path_GBM(const OptionSpec& opt, const std::vector<double>& Z, Tape* tape=nullptr) {
        const int n = opt.steps;
        const double dt = opt.T / n;
        const double sqdt = std::sqrt(dt);

        T S, r, sigma;
        if constexpr (std::is_same_v<T, ADouble>) {
            S     = make_var(opt.S0,  tape);
            r     = make_var(opt.r,   tape);
            sigma = make_var(opt.sigma, tape);
        } else {
            S     = T(opt.S0);
            r     = T(opt.r);
            sigma = T(opt.sigma);
        }

        std::vector<T> path; path.reserve(n);
        for (int i = 0; i < n; ++i) {
            T drift = (r - 0.5 * sigma * sigma - opt.q) * dt;
            T diff  = sigma * sqdt * Z[i];
            S = S * exp(drift + diff);
            path.push_back(S);
        }
        return path;
    }

    /** @brief Simulate a GBM path from explicit active/basic variables.
     *  @tparam T double or ADouble
     *  @param opt Option spec
     *  @param Z   Standard normal shocks
     *  @param S0v Active/basic initial spot
     *  @param rv  Active/basic risk-free rate
     *  @param sigv Active/basic volatility
     */
    template <typename T>
    std::vector<T> simulate_path_GBM_from_vars(const OptionSpec& opt, const std::vector<double>& Z, const T& S0v, const T& rv, const T& sigv) {
        const int n = opt.steps; const double dt = opt.T / n; const double sqdt = std::sqrt(dt);
        std::vector<T> path; path.reserve(n);
        T S = S0v;
        for (int i = 0; i < n; ++i) {
            T drift = (rv - 0.5 * sigv * sigv - opt.q) * dt;
            T diff  = sigv * sqdt * Z[i];
            S = S * exp(drift + diff);
            path.push_back(S);
        }
        return path;
    }

} 
#endif
