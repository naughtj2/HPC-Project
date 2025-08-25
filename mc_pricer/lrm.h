/**
 * @file lrm.h
 * @brief Likelihood Ratio Method (LRM) estimators.
 */
#ifndef PRICER_LRM_H
#define PRICER_LRM_H

#include "options.h"
#include "rng.h"
#include <tuple>

namespace pricer {

    /**
 * @brief Likelihood Ratio Method (score function) estimators for greeks.
 */
    struct LRM {
        /** @brief Single-path LRM estimators for delta and vega of a European option.
         *  @param opt Option specification
         *  @param rng Random number generator
         *  @param antithetic Whether to use antithetic variates
         *  @param is_put If true, returns put estimators; otherwise call
         *  @return Pair (delta_estimator, vega_estimator)
         */
        static std::tuple<double,double> delta_vega_european(const OptionSpec& opt, RNG& rng, bool antithetic, bool is_put);
    };

} 

#endif
