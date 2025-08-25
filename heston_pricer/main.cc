/**
 * @file main.cc
 * @brief CLI entry point for Monte Carlo option pricer and greeks.
 */
#include <iostream>
#include <iomanip>
#include <chrono>
#include "args.h"
#include "bs_analytic.h"
using namespace pricer;
/** 
 * @brief Program entry point.
 * Parses CLI args, runs pricing (multithreaded, optional QMC/control variate),
 * computes finite-difference greeks, and prints a summary.
 * @param argc Argument count
 * @param argv Argument vector
 * @return 0 on success
 */
int main(int argc, char** argv){
    std::ios::sync_with_stdio(false);
    std::cout.setf(std::ios::fixed); std::cout<<std::setprecision(6);
    Args args; args.parse(argc, argv);
    
    auto t0=std::chrono::high_resolution_clock::now();
    const unsigned long long N=args.cfg.paths; const int Tn=std::max(1,args.cfg.threads); const unsigned long long chunk=(N+Tn-1)/Tn;
    
    std::vector<std::thread> threads; std::vector<ThreadSums> partial(Tn);
    PathFunc path_func = make_path_functor(args.opt, args.model, args.H);
    
    for(int t=0;t<Tn;++t){
        unsigned long long i0=std::min<unsigned long long>(N, t*chunk), i1=std::min<unsigned long long>(N, (unsigned long long)(t+1)*chunk);
        if(i0>=i1){ partial[t]=ThreadSums{}; continue; }
        threads.emplace_back([&,i0,i1,t]{ simulate_batch(i0,i1,args.opt,args.H,args.model,args.cfg,path_func,partial[t], 0x9E3779B97F4A7C15ULL*(t+17)); });
    }

    for(auto& th: threads) th.join();
    
    ThreadSums tot{};
    for(const auto& p: partial){
        tot.sum+=p.sum; tot.sum2+=p.sum2; tot.sumC+=p.sumC; tot.sumC2+=p.sumC2; tot.sumPC+=p.sumPC;
        tot.sum_dS+=p.sum_dS; tot.sum_dS2+=p.sum_dS2; tot.sum_dR+=p.sum_dR; tot.sum_dR2+=p.sum_dR2;
        tot.sum_dSigma+=p.sum_dSigma; tot.sum_dSigma2+=p.sum_dSigma2;
        tot.sum_delta_lrm+=p.sum_delta_lrm; tot.sum_delta_lrm2+=p.sum_delta_lrm2;
        tot.sum_vega_lrm+=p.sum_vega_lrm;   tot.sum_vega_lrm2+=p.sum_vega_lrm2;
    }

    MCResult res = combine_results_and_price(args.opt, args.cfg, tot);
    auto t1=std::chrono::high_resolution_clock::now();
    MCResult fd = finite_differences(args.opt, args.H, args.model, args.cfg);
    res.delta_fd=fd.delta_fd; res.gamma_fd=fd.gamma_fd; res.vega_fd=fd.vega_fd;
    auto t2=std::chrono::high_resolution_clock::now();
    
    std::cout<<"Model: "<<(args.model==Model::BS?"Black–Scholes (GBM)":"Heston (Euler full truncation)")<<"\n";
    std::cout<<"Option: ";
    
    switch(args.opt.type){
        case OptionType::EuropeanCall: std::cout<<"European Call"; break;
        case OptionType::EuropeanPut:  std::cout<<"European Put"; break;
        case OptionType::AsianArithmeticCall: std::cout<<"Asian Arithmetic Call"; break;
        case OptionType::AsianArithmeticPut:  std::cout<<"Asian Arithmetic Put"; break;
        case OptionType::AsianAvgStrikeCall:  std::cout<<"Asian Average-Strike Call"; break;
        case OptionType::AsianAvgStrikePut:   std::cout<<"Asian Average-Strike Put"; break;
    }
    
    std::cout << "\nS0=" << args.opt.S0 << ", K=" << args.opt.K << ", r=" << args.opt.r << ", q=" << args.opt.q
              << ", T=" << args.opt.T << ", steps=" << args.opt.steps << ", paths=" << args.cfg.paths << "\n";
    
    if(args.model==Model::BS) std::cout<<"sigma="<<args.opt.sigma<<"\n";
    
    else std::cout<<"kappa="<<args.H.kappa<<", theta="<<args.H.theta<<", xi="<<args.H.xi<<", v0="<<args.H.v0<<", rho="<<args.H.rho<<"\n";
    
    std::cout<<"antithetic="<<(args.cfg.antithetic?"on":"off")<<", qmc(Halton)="<<(args.cfg.use_qmc?"on":"off")
    <<", control_variate(S_T)="<<(args.cfg.use_control_variate?"on":"off")<<", threads="<<args.cfg.threads<<"\n\n";
    
    std::cout<<"Price  = "<<res.price<<"  (SE "<<res.stderr<<")\n\n";
    
    std::cout<<"Greeks:  Finite Differences (CRN):\n"<<"  Delta = "<<res.delta_fd<<"\n"<<"  Gamma = "<<res.gamma_fd<<"\n"
    <<"  "<<(args.model==Model::BS?"Vega ":"dPrice/dv0")<<" = "<<res.vega_fd<<"\n";
    
    if(args.model==Model::BS && (args.opt.type==OptionType::EuropeanCall || args.opt.type==OptionType::EuropeanPut)){
        double priceA = (args.opt.type==OptionType::EuropeanCall) ? BSAnalytic::call(args.opt.S0,args.opt.K,args.opt.r,args.opt.q,args.opt.sigma,args.opt.T)
                                                                  : BSAnalytic::put (args.opt.S0,args.opt.K,args.opt.r,args.opt.q,args.opt.sigma,args.opt.T);
        
        double deltaA = (args.opt.type==OptionType::EuropeanCall) ? BSAnalytic::delta_call(args.opt.S0,args.opt.K,args.opt.r,args.opt.q,args.opt.sigma,args.opt.T)
                                                                  : BSAnalytic::delta_put (args.opt.S0,args.opt.K,args.opt.r,args.opt.q,args.opt.sigma,args.opt.T);
        double vegaA  = BSAnalytic::vega(args.opt.S0,args.opt.K,args.opt.r,args.opt.q,args.opt.sigma,args.opt.T);
        std::cout<<"\nGreeks : Likelihood Ratio Method (BS European):\n"
                 <<"  Delta = "<< (tot.sum_delta_lrm/(double)args.cfg.paths) <<"\n"
                 <<"  Vega  = "<< (tot.sum_vega_lrm /(double)args.cfg.paths) <<"\n";
        std::cout<<"\n(Analytic BS check: price="<<priceA<<", delta="<<deltaA<<", vega="<<vegaA<<")\n";
    }
    
    res.delta_aad = tot.sum_dS/(double)args.cfg.paths; 
    res.rho_aad = tot.sum_dR/(double)args.cfg.paths; 
    res.vega_aad = tot.sum_dSigma/(double)args.cfg.paths;

    std::cout<< "\n  Greeks: AAD (pathwise):\n" <<"  Delta = "<<res.delta_aad<<"\n" <<"  "<<(args.model==Model::BS?"Vega ":"dPrice/dv0")<<" = "<<res.vega_aad<<"\n" <<"  Rho   = "<<res.rho_aad<<"\n";
    
    auto ms1=std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count(); auto ms2=std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    
    std::cout<<"\nTiming: price+AAD(+LRM) pass = "<<ms1<<" ms,  finite-diff pass = "<<ms2<<" ms\n";
    return 0;
}
