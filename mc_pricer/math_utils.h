/**
 * @file math_utils.h
 * @brief Math utilities: Normal PDF/CDF and inverse CDF.
 */
#ifndef PRICER_MATH_UTILS_H
#define PRICER_MATH_UTILS_H

#include <cmath>
#include <limits>

namespace pricer {

    /** @brief Pi constant. */
    constexpr double PI = 3.141592653589793238462643383279502884;

    /** @brief Standard normal CDF \f$\Phi(x)\f$. */
    inline double phi_cdf(double x) {
        return 0.5 * std::erfc(-x / std::sqrt(2.0));
    }

    /** @brief Standard normal PDF \f$\varphi(x)\f$. */
    inline double phi_pdf(double x) {
        return std::exp(-0.5 * x * x) / std::sqrt(2.0 * PI);
    }

    /** @brief Approximate inverse standard normal CDF.
     * @param p Probability in (0,1)
     * @return x such that Phi(x) \approx p
     */
    inline double inv_phi(double p) {
        if (p <= 0.0) return -std::numeric_limits<double>::infinity();
        if (p >= 1.0) return  std::numeric_limits<double>::infinity();

        static const double a1 = -3.969683028665376e+01;
        static const double a2 =  2.209460984245205e+02;
        static const double a3 = -2.759285104469687e+02;
        static const double a4 =  1.383577518672690e+02;
        static const double a5 = -3.066479806614716e+01;
        static const double a6 =  2.506628277459239e+00;

        static const double b1 = -5.447609879822406e+01;
        static const double b2 =  1.615858368580409e+02;
        static const double b3 = -1.556989798598866e+02;
        static const double b4 =  6.680131188771972e+01;
        static const double b5 = -1.328068155288572e+01;

        static const double c1 = -7.784894002430293e-03;
        static const double c2 = -3.223964580411365e-01;
        static const double c3 = -2.400758277161838e+00;
        static const double c4 = -2.549732539343734e+00;
        static const double c5 =  4.374664141464968e+00;
        static const double c6 =  2.938163982698783e+00;

        static const double d1 =  7.784695709041462e-03;
        static const double d2 =  3.224671290700398e-01;
        static const double d3 =  2.445134137142996e+00;
        static const double d4 =  3.754408661907416e+00;

        static const double plow  = 0.02425;
        static const double phigh = 1 - plow;

        double q, r, x;
        if (p < plow) {
            q = std::sqrt(-2 * std::log(p));
            x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
                ((((d1*q + d2)*q + d3)*q + d4)*q + 1);
        } else if (p > phigh) {
            q = std::sqrt(-2 * std::log(1 - p));
            x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
                ((((d1*q + d2)*q + d3)*q + d4)*q + 1);
        } else {
            q = p - 0.5;
            r = q * q;
            x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6) * q /
                (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1);
        }

        double e = phi_cdf(x) - p;
        double u = e / (std::sqrt(2.0 * PI) * std::exp(-0.5 * x * x));
        x -= u / (1 + 0.5 * x * u);
        return x;
    }

} 

#endif
