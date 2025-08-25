/**
 * @file options.h
 * @brief Option types, model selection, and payoff utilities.
 */
#ifndef PRICER_OPTIONS_H
#define PRICER_OPTIONS_H

#include <vector>
#include "aad.h"

namespace pricer {
    /** @brief Supported option contracts. */
    enum class OptionType { EuropeanCall, EuropeanPut, AsianArithmeticCall, AsianArithmeticPut, AsianAvgStrikeCall, AsianAvgStrikePut };
    /** @brief Supported models. */
    enum class Model { BS, Heston };
    struct OptionSpec { 
        OptionType type{OptionType::AsianArithmeticCall}; 
        double S0{100.0}, K{100.0}, r{0.02}, q{0.0}, sigma{0.20}, T{1.0}; int steps{252}; 
    };
    /**
 * @brief Heston model parameters.
 */
    struct HestonSpec { double kappa{5.0}, theta{0.05}, xi{0.5}, v0{0.05}, rho{-0.8}; };
    template <typename T>
    /** @brief Generic payoff for the supported contracts.
  *  @tparam T double or ADouble
  *  @param opt Option spec
  *  @param path Simulated path (last is S_T)
  *  @return Present value of payoff
  */
    T payoff(const OptionSpec& opt, const std::vector<T>& path){
        const T& ST=path.back(); double disc=std::exp(-opt.r*opt.T); auto zeroT=asT(0.0,ST);
        if(opt.type==OptionType::EuropeanCall) return disc * fmax_generic(ST - asT(opt.K,ST), zeroT);
        if(opt.type==OptionType::EuropeanPut ) return disc * fmax_generic(asT(opt.K,ST) - ST, zeroT);
        const T& ref=path.front(); T avg=asT(0.0,ref); for(const auto& s: path) avg=avg+s; avg=avg/(double)path.size();
        if(opt.type==OptionType::AsianArithmeticCall) return disc * fmax_generic(avg - asT(opt.K,ref), zeroT);
        if(opt.type==OptionType::AsianArithmeticPut ) return disc * fmax_generic(asT(opt.K,ref) - avg, zeroT);
        if(opt.type==OptionType::AsianAvgStrikeCall)  return disc * fmax_generic(ST - avg, zeroT);
        if(opt.type==OptionType::AsianAvgStrikePut )  return disc * fmax_generic(avg - ST, zeroT);
        return asT(0.0, ST);
    }
} 
#endif
