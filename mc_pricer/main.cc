/**
 * @file main.cc
 * @brief CLI entry point for Monte Carlo option pricer.
 */
#include <iostream>
#include <iomanip>
#include <chrono>
#include "args.h"
#include "bs_analytic.h"

using namespace pricer;

/**
 * @brief Program entry point.
 *
 * Parses CLI arguments, runs Monte Carlo pricing, finite-difference greeks,
 * and (when applicable) AAD and LRM estimators, then prints a report.
 *
 * @param argc Argument count
 * @param argv Argument vector
 * @return Exit code 
 */
int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);

    Args args; args.parse(argc, argv);

    auto t0 = std::chrono::high_resolution_clock::now();

    MCResult res = monte_carlo_price(args.opt, args.cfg);

    auto t1 = std::chrono::high_resolution_clock::now();

    MCResult fd = finite_differences(args.opt, args.cfg, res.price);
    res.delta_fd = fd.delta_fd; res.gamma_fd = fd.gamma_fd; res.vega_fd = fd.vega_fd;

    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Option: ";
    switch (args.opt.type) {
        case OptionType::EuropeanCall: std::cout << "European Call"; break;
        case OptionType::EuropeanPut: std::cout << "European Put"; break;
        case OptionType::AsianArithmeticCall: std::cout << "Asian Arithmetic Call"; break;
        case OptionType::AsianArithmeticPut: std::cout << "Asian Arithmetic Put"; break;
        case OptionType::AsianAvgStrikeCall: std::cout << "Asian Average-Strike Call"; break;
        case OptionType::AsianAvgStrikePut: std::cout << "Asian Average-Strike Put"; break;
    }
    std::cout << "\nS0=" << args.opt.S0 << ", K=" << args.opt.K << ", r=" << args.opt.r
              << ", sigma=" << args.opt.sigma << ", T=" << args.opt.T
              << ", steps=" << args.opt.steps << ", paths=" << args.cfg.paths << "\n";
    std::cout << "antithetic=" << (args.cfg.antithetic?"on":"off")
              << ", qmc(Halton)=" << (args.cfg.use_qmc?"on":"off")
              << ", control_variate(S_T)=" << (args.cfg.use_control_variate?"on":"off")
              << ", threads=" << args.cfg.threads << "\n\n";

    std::cout << "Price  = " << res.price << "  (SE " << res.stderr << ")\n";

    std::cout << "\nGreeks:  Finite Differences (CRN):\n";
    std::cout << "  Delta = " << res.delta_fd << "\n";
    std::cout << "  Gamma = " << res.gamma_fd << "\n";
    std::cout << "  Vega  = " << res.vega_fd << "\n";

    if (args.opt.type == OptionType::EuropeanCall || args.opt.type == OptionType::EuropeanPut) {
        std::cout << "\nGreeks : Likelihood Ratio Method (European):\n";
        std::cout << "  Delta = " << res.delta_lrm << "\n";
        std::cout << "  Vega  = " << res.vega_lrm << "\n";

        double priceA = (args.opt.type == OptionType::EuropeanCall)
            ? BSAnalytic::call(args.opt.S0, args.opt.K, args.opt.r, args.opt.q, args.opt.sigma, args.opt.T)
            : BSAnalytic::put (args.opt.S0, args.opt.K, args.opt.r, args.opt.q, args.opt.sigma, args.opt.T);
        double deltaA = (args.opt.type == OptionType::EuropeanCall)
            ? BSAnalytic::delta_call(args.opt.S0, args.opt.K, args.opt.r, args.opt.q, args.opt.sigma, args.opt.T)
            : BSAnalytic::delta_put (args.opt.S0, args.opt.K, args.opt.r, args.opt.q, args.opt.sigma, args.opt.T);
        double vegaA = BSAnalytic::vega(args.opt.S0, args.opt.K, args.opt.r, args.opt.q, args.opt.sigma, args.opt.T);

        std::cout << "\n(Analytic BS check: price=" << priceA
                  << ", delta=" << deltaA
                  << ", vega="  << vegaA  << ")\n";
    }

    std::cout << "\nGreeks: AAD (pathwise):\n";
    std::cout << "  Delta = " << res.delta_aad << "\n";
    std::cout << "  Vega  = " << res.vega_aad  << "\n";
    std::cout << "  Rho   = " << res.rho_aad   << "\n";

    auto ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    auto ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "\nTiming: price+AAD+LRM pass = " << ms1
              << " ms,  finite-diff pass = " << ms2 << " ms\n";

    return 0;
}
