/**
 * @file bs_analytic.h
 * @brief Black–Scholes formulas and greeks.
 */
#ifndef PRICER_BS_ANALYTIC_H
#define PRICER_BS_ANALYTIC_H

namespace pricer {
    /**
     * @brief Static utilities for Black–Scholes pricing and greeks.
     */
    struct BSAnalytic {
        /** @brief d1 term of Black–Scholes. */
        static double d1(double S,double K,double r,double q,double sigma,double T);
        /** @brief European call price. */
        static double call(double S,double K,double r,double q,double sigma,double T);
        /** @brief European put price. */
        static double put (double S,double K,double r,double q,double sigma,double T);
        /** @brief Delta of a European call. */
        static double delta_call(double S,double K,double r,double q,double sigma,double T);
        /** @brief Delta of a European put. */
        static double delta_put (double S,double K,double r,double q,double sigma,double T);
        /** @brief Vega of a European option. */
        static double vega(double S,double K,double r,double q,double sigma,double T);
    };
} 
#endif
