/**
 * @file options.h
 * @brief Option types and payoffs used by the pricer.
 */
#ifndef PRICER_OPTIONS_H
#define PRICER_OPTIONS_H

#include <vector>
#include "aad.h"

namespace pricer {

    /** @brief Supported option contracts. */
    enum class OptionType { EuropeanCall, EuropeanPut, AsianArithmeticCall, AsianArithmeticPut, AsianAvgStrikeCall, AsianAvgStrikePut };

    /**
 * @brief Specification of an option and market/contract parameters.
 * @var S0 Spot price at t=0
 * @var K  Strike
 * @var r  Risk-free continuously compounded rate
 * @var q  Continuous dividend yield
 * @var sigma Volatility
 * @var T  Time to maturity in years
 * @var steps Time steps along the path
 */
    struct OptionSpec {
        OptionType type{OptionType::AsianArithmeticCall};
        double S0{100.0}, K{100.0}, r{0.02}, q{0.0}, sigma{0.20}, T{1.0};
        int steps{252};
    };

    /** @brief Generic payoff function for all supported options.
     *  @tparam T numeric or AAD type exposing +,-,*,/,exp,fmax_generic
     *  @param opt Option specification
     *  @param path Simulated price path (last element is S_T)
     *  @return Present value of payoff
     */
    template <typename T>
    T payoff(const OptionSpec& opt, const std::vector<T>& path) {
        const T& ST = path.back();
        double disc = std::exp(-opt.r * opt.T);
        auto zeroT = asT(0.0, ST);

        if (opt.type == OptionType::EuropeanCall) {
            return disc * fmax_generic(ST - asT(opt.K, ST), zeroT);
        } else if (opt.type == OptionType::EuropeanPut) {
            return disc * fmax_generic(asT(opt.K, ST) - ST, zeroT);
        } else if (opt.type == OptionType::AsianArithmeticCall || opt.type == OptionType::AsianArithmeticPut) {
            const T& ref = path.front();
            T avg = asT(0.0, ref);
            for (const auto& s : path) avg = avg + s;
            avg = avg / (double)path.size();
            if (opt.type == OptionType::AsianArithmeticCall)
                return disc * fmax_generic(avg - asT(opt.K, ref), zeroT);
            else
                return disc * fmax_generic(asT(opt.K, ref) - avg, zeroT);
        } else if (opt.type == OptionType::AsianAvgStrikeCall || opt.type == OptionType::AsianAvgStrikePut) {
            const T& ref = path.front();
            T avg = asT(0.0, ref);
            for (const auto& s : path) avg = avg + s;
            avg = avg / (double)path.size();
            if (opt.type == OptionType::AsianAvgStrikeCall)
                return disc * fmax_generic(ST - avg, zeroT);
            else
                return disc * fmax_generic(avg - ST, zeroT);
        }
        return asT(0.0, ST);
    }

} 

#endif
