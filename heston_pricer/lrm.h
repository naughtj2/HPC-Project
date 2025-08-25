/**
 * @file lrm.h
 * @brief Likelihood Ratio Method (score function) estimators.
 */
#ifndef PRICER_LRM_H
#define PRICER_LRM_H

#include "options.h"
#include "rng.h"
#include <tuple>
#include <algorithm>
#include <cmath>

namespace pricer {
    /**
     * @brief LRM estimators for greeks under Black–Scholes.
     */
    struct LRM {
        /** @brief Single-path LRM estimators for delta and vega of a European option.
    *  @return Pair (delta_est, vega_est)
    */
        static std::tuple<double,double> delta_vega_european_bs(const OptionSpec& opt, RNG& rng, bool antithetic);
    };
} 
#endif
