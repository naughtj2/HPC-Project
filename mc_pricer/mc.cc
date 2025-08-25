/**
 * @file mc.cc
 * @brief Monte Carlo driver and finite-difference pass implementations.
 */
#include "mc.h"

namespace pricer {

    /** @brief Monte Carlo price with optional control variate, QMC, antithetic, multithreading. */
    MCResult monte_carlo_price(const OptionSpec& opt, const MCConfig& cfg) {
        const unsigned long long N = cfg.paths;
        const int Tn = std::max(1, cfg.threads);
        const unsigned long long chunk = (N + Tn - 1) / Tn;

        std::vector<std::thread> threads;
        std::vector<ThreadSums> partial(Tn);

        auto path_func = make_path_functor(opt);

        for (int t = 0; t < Tn; ++t) {
            unsigned long long i0 = std::min<unsigned long long>(N, t * chunk);
            unsigned long long i1 = std::min<unsigned long long>(N, (unsigned long long)(t+1) * chunk);
            if (i0 >= i1) { partial[t] = ThreadSums{}; continue; }
            threads.emplace_back([&, i0, i1, t]{
                simulate_batch(i0, i1, opt, cfg, path_func, partial[t], 0x9E3779B97F4A7C15ULL * t);
            });
        }
        for (auto& th : threads) th.join();

        ThreadSums tot{};
        for (const auto& p : partial) {
            tot.sum += p.sum; tot.sum2 += p.sum2;
            tot.sumC += p.sumC; tot.sumC2 += p.sumC2; tot.sumPC += p.sumPC;
            tot.sum_dS += p.sum_dS; tot.sum_dS2 += p.sum_dS2;
            tot.sum_dR += p.sum_dR; tot.sum_dR2 += p.sum_dR2;
            tot.sum_dSigma += p.sum_dSigma; tot.sum_dSigma2 += p.sum_dSigma2;
            tot.sum_delta_lrm += p.sum_delta_lrm; tot.sum_delta_lrm2 += p.sum_delta_lrm2;
            tot.sum_vega_lrm  += p.sum_vega_lrm;  tot.sum_vega_lrm2  += p.sum_vega_lrm2;
        }

        double price_raw = tot.sum / (double)N;
        double var_raw = std::max(0.0, tot.sum2 / (double)N - price_raw * price_raw);

        double price = price_raw;
        double variance = var_raw;

        if (cfg.use_control_variate) {
            double C_bar = tot.sumC / (double)N;
            double varC  = std::max(0.0, tot.sumC2 / (double)N - C_bar * C_bar);
            double covPC = tot.sumPC / (double)N - price_raw * C_bar;
            double beta = (varC > 0 ? covPC / varC : 0.0);
            double EC = opt.S0 * std::exp(-opt.q * opt.T);
            price = price_raw - beta * (C_bar - EC);
            variance = var_raw + beta*beta*varC - 2*beta*covPC;
        }

        MCResult res; res.price = price; res.stderr = std::sqrt(variance / (double)N);

        attach_aad_and_lrm(opt, cfg, tot, res);
        return res;
    }

    /** @brief Second pass to estimate delta, gamma, vega with CRN finite differences. */
    MCResult finite_differences(const OptionSpec& opt, const MCConfig& cfg, double /*base_price*/) {
        const unsigned long long N = cfg.paths;
        const int Tn = std::max(1, cfg.threads);
        const unsigned long long chunk = (N + Tn - 1) / Tn;

        std::vector<std::thread> threads;
        std::mutex mtx;

        struct Acc { double dS{0}, dS2{0}, g{0}, v{0}; } acc;

        for (int t = 0; t < Tn; ++t) {
            unsigned long long i0 = std::min<unsigned long long>(N, t * chunk);
            unsigned long long i1 = std::min<unsigned long long>(N, (unsigned long long)(t+1) * chunk);
            if (i0 >= i1) continue;
            threads.emplace_back([&, i0, i1]{
                RNG rng(cfg.seed + 0xABCDEF1234567890ULL + i0*1315423911ULL);
                const int n = opt.steps; const double dt = opt.T / n; const double sqdt = std::sqrt(dt);
                std::vector<double> Z(n), u(n);
                Halton halton(n);

                auto gen_normals = [&](unsigned long long path_index){
                    if (cfg.use_qmc) {
                        static thread_local std::vector<double> shift;
                        if (shift.size() != (size_t)n) {
                            shift.resize(n);
                            for (int k = 0; k < n; ++k) shift[k] = rng.uniform();
                        }
                        halton.point((unsigned long long)path_index + 1ULL, u);
                        for (int k = 0; k < n; ++k) { double uu = u[k] + shift[k]; uu -= std::floor(uu); Z[k] = inv_phi(uu); }
                    } else {
                        for (int k = 0; k < n; ++k) Z[k] = rng.normal();
                    }
                };

                auto makeF = [&](const OptionSpec& o){ return make_path_functor(o); };

                OptionSpec upS=opt, dnS=opt, upV=opt, dnV=opt;
                double hS = opt.S0 * cfg.bump_rel_S;
                double hV = opt.sigma * cfg.bump_rel_sigma;
                upS.S0 += hS; dnS.S0 -= hS;
                upV.sigma += hV; dnV.sigma -= hV; if (dnV.sigma < 1e-8) dnV.sigma = 1e-8;

                auto fopt = makeF(opt);
                auto fupS = makeF(upS);
                auto fdnS = makeF(dnS);
                auto fupV = makeF(upV);
                auto fdnV = makeF(dnV);

                double sum_dS=0, sum_dS2=0, sum_g=0, sum_v=0;
                for (unsigned long long i = i0; i < i1; ++i) {
                    gen_normals(i);
                    auto one = [&](auto& F){
                        double pv = F(Z).first;
                        if (cfg.antithetic) {
                            for (double& z:Z) z=-z;
                            pv = 0.5*(pv + F(Z).first);
                        }
                        return pv;
                    };
                    double p_base = one(fopt);
                    double p_upS = one(fupS);
                    double p_dnS = one(fdnS);
                    double p_upV = one(fupV);
                    double p_dnV = one(fdnV);

                    double dS = (p_upS - p_dnS) / (2*hS);
                    double g  = (p_upS - 2*p_base + p_dnS) / (hS*hS);
                    double v  = (p_upV - p_dnV) / (2*hV);

                    sum_dS += dS; sum_dS2 += dS*dS;
                    sum_g  += g;  sum_v  += v;
                }

                std::scoped_lock lk(mtx);
                acc.dS  += sum_dS; acc.dS2 += sum_dS2; acc.g += sum_g; acc.v += sum_v;
            });
        }
        for (auto& th : threads) th.join();

        MCResult out; out.delta_fd = acc.dS / (double)cfg.paths; out.gamma_fd = acc.g / (double)cfg.paths; out.vega_fd = acc.v / (double)cfg.paths;
        return out;
    }

    /** @brief Extract AAD pathwise greeks and LRM (for European) from thread aggregates. */
    void attach_aad_and_lrm(const OptionSpec& opt, const MCConfig& cfg, const ThreadSums& tot, MCResult& res) {
        double N = (double)cfg.paths;
        res.delta_aad = tot.sum_dS / N;
        res.vega_aad  = tot.sum_dSigma / N;
        res.rho_aad   = tot.sum_dR / N;

        if (opt.type == OptionType::EuropeanCall || opt.type == OptionType::EuropeanPut) {
            res.delta_lrm = tot.sum_delta_lrm / N;
            res.vega_lrm  = tot.sum_vega_lrm / N;
        }
    }

} 
