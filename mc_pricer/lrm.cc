/**
 * @file lrm.cc
 * @brief Implementation of LRM estimators for European options.
 */
#include "lrm.h"
#include <algorithm>
#include <cmath>

namespace pricer {

    /** @brief One-sample LRM estimators for delta and vega (European).
     */
    std::tuple<double,double>
    LRM::delta_vega_european(const OptionSpec& opt, RNG& rng, bool antithetic, bool is_put)
    {
        double T = opt.T, sigma = opt.sigma, r = opt.r, q = opt.q, S0 = opt.S0, K = opt.K;
        double vsqrtT = sigma * std::sqrt(T);

        auto one = [&](double Z){
            double ST = S0 * std::exp((r - q - 0.5*sigma*sigma)*T + vsqrtT * Z);
            double disc = std::exp(-r*T);
            double payoff = is_put ? std::max(K - ST, 0.0) : std::max(ST - K, 0.0);
            double w_delta = Z / (S0 * vsqrtT);
            double w_vega  = (Z*Z - 1.0) / sigma;
            return std::make_tuple(disc * payoff * w_delta, disc * payoff * w_vega);
        };

        if (!antithetic) {
            double Z = rng.normal();
            return one(Z);
        } else {
            double Z = rng.normal();
            auto [d1,v1] = one(Z);
            auto [d2,v2] = one(-Z);
            return std::make_tuple(0.5*(d1+d2), 0.5*(v1+v2));
        }
    }

} 
