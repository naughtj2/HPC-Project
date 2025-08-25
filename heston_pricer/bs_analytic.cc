/**
 * @file bs_analytic.cc
 * @brief Implementation of Black–Scholes pricing and greeks.
 */
#include "bs_analytic.h"
#include "math_utils.h"
#include <algorithm>
#include <cmath>

namespace pricer {
    /** @brief Compute d1. */
    double BSAnalytic::d1(double S,double K,double r,double q,double sigma,double T){ 
        double vs=sigma*std::sqrt(T); 
        return (std::log(S/K) + (r - q + 0.5*sigma*sigma)*T)/vs; 
    }

    /** @brief Closed-form call price. */
    double BSAnalytic::call(double S,double K,double r,double q,double sigma,double T){ 
        if(sigma<=0||T<=0) return std::max(S*std::exp(-q*T)-K*std::exp(-r*T),0.0); 
        double d_1=d1(S,K,r,q,sigma,T), d_2=d_1-sigma*std::sqrt(T); 
        return S*std::exp(-q*T)*phi_cdf(d_1) - K*std::exp(-r*T)*phi_cdf(d_2); 
    }

    /** @brief Closed-form put price. */
    double BSAnalytic::put (double S,double K,double r,double q,double sigma,double T){ 
        if(sigma<=0||T<=0) return std::max(K*std::exp(-r*T)-S*std::exp(-q*T),0.0);
         double d_1=d1(S,K,r,q,sigma,T), d_2=d_1-sigma*std::sqrt(T); 
         return K*std::exp(-r*T)*phi_cdf(-d_2) - S*std::exp(-q*T)*phi_cdf(-d_1); 
        }

    /** @brief Call delta. */
    double BSAnalytic::delta_call(double S,double K,double r,double q,double sigma,double T){ 
        double d_1=d1(S,K,r,q,sigma,T); 
        return std::exp(-q*T)*phi_cdf(d_1); 
    }
    /** @brief Put delta. */
    double BSAnalytic::delta_put (double S,double K,double r,double q,double sigma,double T){ 
        double d_1=d1(S,K,r,q,sigma,T); 
        return std::exp(-q*T)*(phi_cdf(d_1)-1.0); 
    }
    /** @brief Option vega. */
    double BSAnalytic::vega(double S,double K,double r,double q,double sigma,double T){ 
        double d_1=d1(S,K,r,q,sigma,T); 
        return S*std::exp(-q*T)*std::sqrt(T)*phi_pdf(d_1); 
    }
} 