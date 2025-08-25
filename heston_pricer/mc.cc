/**
 * @file mc.cc
 * @brief Monte Carlo driver (pricing, control variate) and finite-difference greeks.
 */
#include "mc.h"
namespace pricer {
    MCResult combine_results_and_price(const OptionSpec& opt, const MCConfig& cfg, const ThreadSums& tot){
        MCResult res; double N=(double)cfg.paths;
        double price_raw = tot.sum / N; 
        double var_raw = std::max(0.0, tot.sum2/N - price_raw*price_raw);
        double price=price_raw, variance=var_raw;
        if(cfg.use_control_variate){
            double C_bar=tot.sumC/N; 
            double varC=std::max(0.0, tot.sumC2/N - C_bar*C_bar); 
            double covPC=tot.sumPC/N - price_raw*C_bar;
            double beta=(varC>0 ? covPC/varC : 0.0); 
            double EC=opt.S0*std::exp(-opt.q*opt.T);
            price = price_raw - beta*(C_bar - EC); 
            variance = var_raw + beta*beta*varC - 2*beta*covPC;
        }
        res.price=price; 
        res.stderr=std::sqrt(variance/N); 
        return res;
    }
    
    MCResult finite_differences(const OptionSpec& opt_in, const HestonSpec& H_in, Model model, const MCConfig& cfg){
        const unsigned long long N=cfg.paths; 
        const int Tn=std::max(1,cfg.threads); 
        const unsigned long long chunk=(N+Tn-1)/Tn;
        std::vector<std::thread> threads; 
        std::mutex mtx; 
        struct Acc{ double dS{0}, dS2{0}, g{0}, v{0}; } acc;
        for(int t=0;t<Tn;++t){
            unsigned long long i0=std::min<unsigned long long>(N, t*chunk), i1=std::min<unsigned long long>(N, (unsigned long long)(t+1)*chunk);
            if(i0>=i1) continue;
            threads.emplace_back([&, i0,i1]{
                RNG rng(cfg.seed + 0xABCDEF1234567890ULL + i0*1315423911ULL); 
                const int n=opt_in.steps; 
                const int dim=(model==Model::BS? n : 2*n);

                std::vector<double> Z(dim), u(dim); Halton halton(dim);
                auto gen_normals = [&](unsigned long long path_index){
                    if(cfg.use_qmc){ 
                        static thread_local std::vector<double> shift; 
                        if(shift.size()!=(size_t)dim){ 
                            shift.resize(dim); 
                            for(int k=0;k<dim;++k) shift[k]=rng.uniform(); 
                        }
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

                OptionSpec opt=opt_in; 
                HestonSpec H=H_in;
                auto Fopt = make_path_functor(opt, model, H);
                OptionSpec upS=opt, dnS=opt; 
                double hS=opt.S0*cfg.bump_rel_S; 
                upS.S0+=hS; 
                dnS.S0-=hS;
                auto FupS=make_path_functor(upS, model, H), FdnS=make_path_functor(dnS, model, H);
                double hV; 
                OptionSpec ov=opt; 
                HestonSpec Hv=H;
                
                if(model==Model::BS){ 
                    ov.sigma += opt.sigma*cfg.bump_rel_sigma; 
                    hV=opt.sigma*cfg.bump_rel_sigma; 
                }
                else { 
                    Hv.v0 += H.v0*cfg.bump_rel_sigma; 
                    hV=H.v0*cfg.bump_rel_sigma; 
                }
                auto FupV=make_path_functor(ov, model, Hv);
                if(model==Model::BS){ 
                    ov=opt; 
                    ov.sigma -= opt.sigma*cfg.bump_rel_sigma; 
                }
                else { 
                    Hv=H; 
                    Hv.v0 -= H.v0*cfg.bump_rel_sigma; 
                    if(Hv.v0<1e-10) Hv.v0=1e-10; 
                }
                auto FdnV=make_path_functor(opt, model, (model==Model::BS? H : Hv));
                auto one=[&](PathFunc& F){ 
                    auto pr=F(Z); 
                    double pv=pr.first; 
                    if(cfg.antithetic){ 
                        for(double& z:Z) z=-z; 
                        auto pr2=F(Z); 
                        for(double& z:Z) z=-z; 
                        pv=0.5*(pv+pr2.first);
                    } 
                    return pv; 
                };
                double sum_dS=0,sum_dS2=0,sum_g=0,sum_v=0;
                for(unsigned long long i=i0;i<i1;++i){ 
                    gen_normals(i); 
                    double p_upS=one(FupS), p_dnS=one(FdnS), p_upV=one(FupV), p_dnV=one(FdnV), base=one(Fopt);
                    double dS=(p_upS - p_dnS)/(2*hS); 
                    double g=(p_upS - 2*base + p_dnS)/(hS*hS); 
                    double v=(p_upV - p_dnV)/(2*hV);
                    sum_dS+=dS; 
                    sum_dS2+=dS*dS; 
                    sum_g+=g; 
                    sum_v+=v; }
                std::scoped_lock lk(mtx); 
                acc.dS+=sum_dS; 
                acc.dS2+=sum_dS2; 
                acc.g+=sum_g; 
                acc.v+=sum_v;
            }
        );
        
    }
        for(auto& th: threads) th.join();
        MCResult out; 
        out.delta_fd=acc.dS/(double)cfg.paths; 
        out.gamma_fd=acc.g/(double)cfg.paths; 
        out.vega_fd=acc.v/(double)cfg.paths; 
        return out;
    }
} 