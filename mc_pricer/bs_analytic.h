/**
 * @file bs_analytic.h
 * @brief Closed-form Black–Scholes(-Merton) formulas and greeks.
 */
#ifndef PRICER_BS_ANALYTIC_H
#define PRICER_BS_ANALYTIC_H

namespace pricer {

    /**
 * @brief Static Black–Scholes pricing and greeks for European options.
 */
    struct BSAnalytic {
        /** @brief Black–Scholes d1 term. */
        static double d1(double S, double K, double r, double q, double sigma, double T);
        /** @brief Black–Scholes d2 term. */
        static double d2(double S, double K, double r, double q, double sigma, double T);
        /** @brief European call price under Black–Scholes. */
        static double call(double S, double K, double r, double q, double sigma, double T);
        /** @brief European put price under Black–Scholes. */
        static double put(double S, double K, double r, double q, double sigma, double T);
        /** @brief Delta of a European call. */
        static double delta_call(double S, double K, double r, double q, double sigma, double T);
        /** @brief Delta of a European put. */
        static double delta_put(double S, double K, double r, double q, double sigma, double T);
        /** @brief Vega of a European option. */
        static double vega(double S, double K, double r, double q, double sigma, double T);
    };

} 

#endif
