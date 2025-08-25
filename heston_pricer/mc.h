/**
 * @file mc.h
 * @brief Monte Carlo engine interfaces, config, and helpers.
 */
#ifndef PRICER_MC_H
#define PRICER_MC_H

#include <vector>
#include <thread>
#include <mutex>
#include <functional>
#include <numeric>
#include <optional>
#include <utility>
#include <limits>
#include "rng.h"
#include "halton.h"
#include "math_utils.h"
#include "options.h"
#include "simulate.h"
#include "aad.h"
#include "lrm.h"

namespace pricer {
        /**
     * @brief Aggregate results from Monte Carlo runs (price, stderr, greeks).
     */
struct MCResult { 
        double price{0}, stderr{0}; 
        double delta_fd{std::numeric_limits<double>::quiet_NaN()}; 
        double gamma_fd{std::numeric_limits<double>::quiet_NaN()}; 
        double vega_fd{std::numeric_limits<double>::quiet_NaN()}; 
        double delta_lrm{std::numeric_limits<double>::quiet_NaN()}; 
        double vega_lrm{std::numeric_limits<double>::quiet_NaN()}; 
        double delta_aad{std::numeric_limits<double>::quiet_NaN()}; 
        double vega_aad{std::numeric_limits<double>::quiet_NaN()}; 
        double rho_aad{std::numeric_limits<double>::quiet_NaN()}; 
    };

        /**
     * @brief Configuration for Monte Carlo simulation (paths, RNG, QMC, CV, bumps).
     */
struct MCConfig { 
        unsigned long long paths{200000}; 
        int threads{(int)std::max(1u,std::thread::hardware_concurrency())}; 
        bool antithetic{true}; 
        bool use_qmc{false}; 
        bool use_control_variate{true}; 
        unsigned long long seed{0xC0FFEEULL}; 
        double bump_rel_S{1e-4}; 
        double bump_rel_sigma{1e-4}; 
        double bump_abs_r{1e-4}; 
    };
        /**
     * @brief Per-thread partial sums used for reduction (price, CV, AAD/LRM terms).
     */
struct ThreadSums { 
        double sum{0}, sum2{0}; 
        double sumC{0}, sumC2{0}, sumPC{0}; 
        double sum_dS{0}, sum_dS2{0}; 
        double sum_dSigma{0}, sum_dSigma2{0}; 
        double sum_dR{0}, sum_dR2{0}; 
        double sum_delta_lrm{0}, sum_delta_lrm2{0}; 
        double sum_vega_lrm{0}, sum_vega_lrm2{0}; 
    };

        /** @brief Path payoff/control functor: maps standard normals to (payoff, control).
     */
using PathFunc = std::function<std::pair<double,double>(const std::vector<double>&)>;
    PathFunc make_path_functor(const OptionSpec& opt, Model model, const HestonSpec& H);
    template <typename PathFunctor>
    /** @brief Simulate paths [i0,i1) and accumulate into @p out.
     *  @tparam PathFunctor Callable returning (payoff, control)
     *  @param i0 First path index (inclusive)
     *  @param i1 Last path index (exclusive)
     *  @param opt Option specification
     *  @param H   Heston parameters (ignored for BS)
     *  @param model Model selector (BS/Heston)
     *  @param cfg MC configuration
     *  @param path_func Path functor (payoff/control)
     *  @param out Thread-local accumulators (output)
     *  @param rng_seed_offset Seed offset for independent streams
     */
    void simulate_batch(unsigned long long i0, unsigned long long i1, const OptionSpec& opt, const HestonSpec& H, Model model, 
                        const MCConfig& cfg, PathFunctor&& path_func, ThreadSums& out, unsigned long long rng_seed_offset);
    
    /** @brief Combine thread totals, apply control variate (if enabled), compute price & stderr. */
    MCResult combine_results_and_price(const OptionSpec& opt, const MCConfig& cfg, const ThreadSums& tot);
    /** @brief Estimate greeks via common-random-number finite differences. */
    MCResult finite_differences(const OptionSpec& opt_in, const HestonSpec& H_in, Model model, const MCConfig& cfg);

    inline PathFunc make_path_functor(const OptionSpec& opt, Model model, const HestonSpec& H){
        const int n=opt.steps; 
        const double dt=opt.T/n; 
        const double sqdt=std::sqrt(dt);

        if(model==Model::BS){
            return PathFunc([=](const std::vector<double>& Z){
                double S=opt.S0;
                for(int i=0;i<n;++i){ 
                    double drift=(opt.r - opt.q - 0.5*opt.sigma*opt.sigma)*dt; 
                    double diff=opt.sigma*sqdt*Z[i]; 
                    S*=std::exp(drift+diff); 
                }

                double disc=std::exp(-opt.r*opt.T); 
                double control=disc*S;

                std::vector<double> path; 
                path.reserve(n); 
                S=opt.S0;

                for(int i=0;i<n;++i){ 
                    double drift=(opt.r - opt.q - 0.5*opt.sigma*opt.sigma)*dt; 
                    double diff=opt.sigma*sqdt*Z[i]; 
                    S*=std::exp(drift+diff); 
                    path.push_back(S); 
                }
                
                double pv=0.0; 
                switch(opt.type){

                    case OptionType::EuropeanCall: pv=disc*std::max(path.back()-opt.K,0.0); break;
                    case OptionType::EuropeanPut:  pv=disc*std::max(opt.K-path.back(),0.0); break;
                    case OptionType::AsianArithmeticCall: { double avg=std::accumulate(path.begin(),path.end(),0.0)/path.size(); pv=disc*std::max(avg-opt.K,0.0); break;}
                    case OptionType::AsianArithmeticPut:  { double avg=std::accumulate(path.begin(),path.end(),0.0)/path.size(); pv=disc*std::max(opt.K-avg,0.0); break;}
                    case OptionType::AsianAvgStrikeCall:  { double avg=std::accumulate(path.begin(),path.end(),0.0)/path.size(); pv=disc*std::max(path.back()-avg,0.0); break;}
                    case OptionType::AsianAvgStrikePut:   { double avg=std::accumulate(path.begin(),path.end(),0.0)/path.size(); pv=disc*std::max(avg-path.back(),0.0); break;}
                } 
                return std::make_pair(pv, control);

            }
        );
        } 
        else {
            return PathFunc([=](const std::vector<double>& Z){
                double S=opt.S0, v=H.v0;
                for(int i=0;i<n;++i){ 
                    double W1=Z[2*i], W2=Z[2*i+1]; 
                    double Z1=W1; 
                    double Z2=H.rho*W1 + std::sqrt(std::max(0.0,1.0-H.rho*H.rho))*W2; 
                    double v_pos=std::max(0.0,v);

                    double dS=(opt.r - opt.q)*S*dt + std::sqrt(v_pos)*S*Z1*sqdt; 
                    double dv=H.kappa*(H.theta - v)*dt + H.xi*std::sqrt(v_pos)*Z2*sqdt; S+=dS; 
                    v=std::max(0.0, v+dv); 
                }

                double disc=std::exp(-opt.r*opt.T); 
                double control=disc*S;
                std::vector<double> path; 
                path.reserve(n); 
                S=opt.S0; 
                v=H.v0;

                for(int i=0;i<n;++i){ 
                    double W1=Z[2*i], W2=Z[2*i+1]; 
                    double Z1=W1; 
                    double Z2=H.rho*W1 + std::sqrt(std::max(0.0,1.0-H.rho*H.rho))*W2; 
                    double v_pos=std::max(0.0,v);
                    double dS=(opt.r - opt.q)*S*dt + std::sqrt(v_pos)*S*Z1*sqdt; 
                    double dv=H.kappa*(H.theta - v)*dt + H.xi*std::sqrt(v_pos)*Z2*sqdt; 
                    S+=dS; 
                    v=std::max(0.0, v+dv); 
                    path.push_back(S); 
                }
                double pv=0.0; switch(opt.type){
                    case OptionType::EuropeanCall: pv=disc*std::max(path.back()-opt.K,0.0); break;
                    case OptionType::EuropeanPut:  pv=disc*std::max(opt.K-path.back(),0.0); break;
                    case OptionType::AsianArithmeticCall: { double avg=std::accumulate(path.begin(),path.end(),0.0)/path.size(); pv=disc*std::max(avg-opt.K,0.0); break;}
                    case OptionType::AsianArithmeticPut:  { double avg=std::accumulate(path.begin(),path.end(),0.0)/path.size(); pv=disc*std::max(opt.K-avg,0.0); break;}
                    case OptionType::AsianAvgStrikeCall:  { double avg=std::accumulate(path.begin(),path.end(),0.0)/path.size(); pv=disc*std::max(path.back()-avg,0.0); break;}
                    case OptionType::AsianAvgStrikePut:   { double avg=std::accumulate(path.begin(),path.end(),0.0)/path.size(); pv=disc*std::max(avg-path.back(),0.0); break;}
                } return std::make_pair(pv, control);
            }
        );
        }
    }
    template <typename PathFunctor>
    /** @brief Definition: simulate a batch and fill @p out. */
    void simulate_batch(unsigned long long i0, unsigned long long i1, const OptionSpec& opt, const HestonSpec& H, Model model, const MCConfig& cfg, PathFunctor&& path_func, ThreadSums& out, unsigned long long rng_seed_offset){
        RNG rng(cfg.seed + rng_seed_offset + i0*1315423911ULL); 
        const int n=opt.steps; 
        const int dim=(model==Model::BS? n : 2*n);
        std::vector<double> Z(dim), u(dim); Halton halton(dim);
        auto gen_normals = [&](unsigned long long path_index){
            if(cfg.use_qmc){ static thread_local std::vector<double> shift; 
                if(shift.size()!=(size_t)dim){ 
                    shift.resize(dim); 
                    for(int k=0; k<dim;++k) shift[k]=rng.uniform(); }
                halton.point((unsigned long long)path_index + 1ULL, u);
                for(int k=0;k<dim;++k){ 
                    double uu=u[k]+shift[k]; 
                    uu-=std::floor(uu); 
                    Z[k]=inv_phi(uu); 
                }
            } 
            else { 
                for(int k=0;k<dim;++k) Z[k]=rng.normal(); 
            }
        };
        ThreadSums local;
        for(unsigned long long i=i0;i<i1;++i){
            gen_normals(i);
            auto pr1=path_func(Z); 
            double payoff_val=pr1.first, control=pr1.second;
            if(cfg.antithetic){ 
                for(double& z:Z) z=-z; 
                auto pr2=path_func(Z); 
                for(double& z:Z) z=-z; 
                payoff_val=0.5*(pr1.first+pr2.first); 
                control=0.5*(pr1.second+pr2.second); 
            }
            
            local.sum+=payoff_val; 
            local.sum2+=payoff_val*payoff_val;
            if(cfg.use_control_variate){ 
                local.sumC+=control; 
                local.sumC2+=control*control; 
                local.sumPC+=payoff_val*control; 
            }
            auto aad_for = [&](const std::vector<double>& Zvec){
                Tape tape; tape.clear();
                if(model==Model::BS){
                    ADouble S0=make_var(opt.S0,&tape), r=make_var(opt.r,&tape), sg=make_var(opt.sigma,&tape);
                    auto pathA=simulate_path_GBM_from_vars<ADouble>(opt, Zvec, S0, r, sg);
                    ADouble pv=payoff<ADouble>(opt, pathA); 
                    int w=pv.idx; reverse_ad(tape,w);
                    return std::array<double,4>{ pv.val(), tape.adj[S0.idx], tape.adj[r.idx], tape.adj[sg.idx] };
                } 
                else {
                    ADouble S0=make_var(opt.S0,&tape), r=make_var(opt.r,&tape);
                    ADouble k=make_var(H.kappa,&tape), th=make_var(H.theta,&tape), xi=make_var(H.xi,&tape), v0=make_var(H.v0,&tape), rh=make_var(H.rho,&tape);
                    auto pathA=simulate_path_Heston_from_vars<ADouble>(opt, H, Zvec, S0, r, k, th, xi, v0, rh);
                    ADouble pv=payoff<ADouble>(opt, pathA); int w=pv.idx; reverse_ad(tape,w);
                    return std::array<double,4>{ pv.val(), tape.adj[S0.idx], tape.adj[r.idx], tape.adj[v0.idx] };
                }
            };
            auto a1=aad_for(Z);
            if(cfg.antithetic){ for(double& z:Z) z=-z; 
                auto a2=aad_for(Z); 
                for(double& z:Z) z=-z; 
                a1[0]=0.5*(a1[0]+a2[0]); 
                a1[1]=0.5*(a1[1]+a2[1]); 
                a1[2]=0.5*(a1[2]+a2[2]); 
                a1[3]=0.5*(a1[3]+a2[3]); 
            }
            local.sum_dS+=a1[1]; 
            local.sum_dS2+=a1[1]*a1[1]; 
            local.sum_dR+=a1[2]; 
            local.sum_dR2+=a1[2]*a1[2]; 
            local.sum_dSigma+=a1[3]; 
            local.sum_dSigma2+=a1[3]*a1[3];
            if(model==Model::BS && (opt.type==OptionType::EuropeanCall || opt.type==OptionType::EuropeanPut)){
                auto lrm=LRM::delta_vega_european_bs(opt, rng, cfg.antithetic); double d=std::get<0>(lrm), v=std::get<1>(lrm);
                local.sum_delta_lrm+=d; 
                local.sum_delta_lrm2+=d*d; 
                local.sum_vega_lrm+=v; 
                local.sum_vega_lrm2+=v*v;
            }
        }
        out=local;
    }
} 
#endif
