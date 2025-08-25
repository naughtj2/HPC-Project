/**
 * @file args.cc
 * @brief Implementation of CLI parsing.
 */
#include "args.h"
#include <cstdlib>
namespace pricer {
    /** @brief Convert option token to enum. */
    OptionType Args::option_from_string(const std::string& s){
        if(s=="euro_call") return OptionType::EuropeanCall;
        if(s=="euro_put")  return OptionType::EuropeanPut;
        if(s=="asian_arith_call") return OptionType::AsianArithmeticCall;
        if(s=="asian_arith_put")  return OptionType::AsianArithmeticPut;
        if(s=="asian_avgstrike_call") return OptionType::AsianAvgStrikeCall;
        if(s=="asian_avgstrike_put")  return OptionType::AsianAvgStrikePut;
        return OptionType::AsianArithmeticCall;
    }
    /** @brief Parse argv and set model, option spec, and MC config. */
    void Args::parse(int argc, char** argv){
        for(int i=1;i<argc;++i){
            std::string a=argv[i]; auto need=[&](int){ return i+1<argc? argv[++i]:nullptr; };
            if(a=="--model"){ std::string m=need(1); model = (m=="heston"? Model::Heston : Model::BS); }
            else if(a=="--option") opt.type=option_from_string(need(1));
            else if(a=="--S0") opt.S0=std::stod(need(1));
            else if(a=="--K")  opt.K =std::stod(need(1));
            else if(a=="--r")  opt.r =std::stod(need(1));
            else if(a=="--q")  opt.q =std::stod(need(1));
            else if(a=="--sigma") opt.sigma=std::stod(need(1));
            else if(a=="--T")  opt.T =std::stod(need(1));
            else if(a=="--steps") opt.steps=std::stoi(need(1));
            else if(a=="--paths") cfg.paths=(unsigned long long)std::stoll(need(1));
            else if(a=="--threads") cfg.threads=std::stoi(need(1));
            else if(a=="--seed") cfg.seed=(unsigned long long)std::stoull(need(1));
            else if(a=="--antithetic") cfg.antithetic=std::stoi(need(1))!=0;
            else if(a=="--qmc") cfg.use_qmc=std::stoi(need(1))!=0;
            else if(a=="--control_variate") cfg.use_control_variate=std::stoi(need(1))!=0;
            else if(a=="--bump_rel_S") cfg.bump_rel_S=std::stod(need(1));
            else if(a=="--bump_rel_sigma") cfg.bump_rel_sigma=std::stod(need(1));
            else if(a=="--bump_abs_r") cfg.bump_abs_r=std::stod(need(1));
            else if(a=="--kappa") H.kappa=std::stod(need(1));
            else if(a=="--theta") H.theta=std::stod(need(1));
            else if(a=="--xi")    H.xi   =std::stod(need(1));
            else if(a=="--v0")    H.v0   =std::stod(need(1));
            else if(a=="--rho")   H.rho  =std::stod(need(1));
            else if(a=="--help"){
                std::cout << "\nUsage: ./pricerheston --model <bs|heston> [options]\n\n"
                        << "Common: --option <euro_call|euro_put|asian_arith_call|asian_arith_put|asian_avgstrike_call|asian_avgstrike_put>\n"
                        << "        --S0 --K --r --q --T --steps --paths --threads --antithetic 0/1 --qmc 0/1 --seed <u64>\n"
                        << "BS:     --sigma <vol>\n"
                        << "Heston: --kappa --theta --xi (vol-of-vol) --v0 (init var) --rho\n\n";
                std::exit(0);
        }
    }
}
} 