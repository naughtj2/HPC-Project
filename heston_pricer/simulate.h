/**
 * @file simulate.h
 * @brief Path simulation under GBM and Heston (double or AAD).
 */
#ifndef PRICER_SIMULATE_H
#define PRICER_SIMULATE_H

#include <vector>
#include <type_traits>
#include <cmath>
#include "options.h"
#include "aad.h"
#include "math_utils.h"

namespace pricer {
    template <typename T>
    /** @brief Simulate a GBM path; if @p tape is provided and T=ADouble, parameters are active.
  *  @tparam T double or ADouble
  */
    std::vector<T> simulate_path_GBM(const OptionSpec& opt, const std::vector<double>& Z, Tape* tape=nullptr){
        const int n=opt.steps; 
        const double dt=opt.T/n; 
        const double sqdt=std::sqrt(dt);

        T S,r,sigma; 
        if constexpr (std::is_same_v<T,ADouble>){ 
            S=make_var(opt.S0,tape); 
            r=make_var(opt.r,tape); 
            sigma=make_var(opt.sigma,tape);
        } 
        else { 
            S=T(opt.S0); 
            r=T(opt.r); 
            sigma=T(opt.sigma); 
        }

        std::vector<T> path; 
        path.reserve(n); 
        for(int i=0;i<n;++i){ 
            T drift=(r - 0.5*sigma*sigma - opt.q) * dt; 
            T diff = sigma * sqdt * Z[i]; 
            S = S * exp(drift + diff); 
            path.push_back(S);
        } 
        return path; 
    }

    template <typename T>
    /** @brief Simulate a GBM path from explicit variable values. */
    std::vector<T> simulate_path_GBM_from_vars(const OptionSpec& opt, const std::vector<double>& Z, const T& S0v, const T& rv, const T& sigv){
        const int n=opt.steps; 
        const double dt=opt.T/n; 
        const double sqdt=std::sqrt(dt);
        
        std::vector<T> path; 
        path.reserve(n); 
        T S=S0v; 
        for(int i=0;i<n;++i){ 
            T drift=(rv - 0.5*sigv*sigv - opt.q)*dt; 
            T diff = sigv * sqdt * Z[i]; 
            S = S * exp(drift + diff); 
            path.push_back(S);
        } 
        return path; 
    }

    template <typename T>
    /** @brief Simulate a Heston path (Euler full truncation); T may be ADouble. */
    std::vector<T> simulate_path_Heston(const OptionSpec& opt, const HestonSpec& H, const std::vector<double>& Z, Tape* tape=nullptr){
        const int n=opt.steps; 
        const double dt=opt.T/n; 
        const double sqdt=std::sqrt(dt);

        T S,r,kappa,theta,xi,v,rho; 
        if constexpr (std::is_same_v<T,ADouble>){ 
            S=make_var(opt.S0,tape); 
            r=make_var(opt.r,tape); 
            kappa=make_var(H.kappa,tape); 
            theta=make_var(H.theta,tape); 
            xi=make_var(H.xi,tape); 
            v=make_var(H.v0,tape); 
            rho=make_var(H.rho,tape);
        } 
        else { 
            S=T(opt.S0); 
            r=T(opt.r); 
            kappa=T(H.kappa); 
            theta=T(H.theta); 
            xi=T(H.xi); 
            v=T(H.v0); 
            rho=T(H.rho); 
        }

        std::vector<T> path; 
        path.reserve(n); 
        for(int i=0;i<n;++i){ 
            double W1=Z[2*i], W2=Z[2*i+1]; 
            T Z1=asT(W1,S);
            T Z2 = rho*Z1 + sqrt( fmax_generic(asT(1.0,S) - rho*rho, asT(0.0,S)) ) * asT(W2,S);
            T v_pos=fmax_generic(v, asT(0.0,v)); 
            T dS=(r - opt.q)*S*dt + sqrt(v_pos)*S*Z1*sqdt; 
            T dv=kappa*(theta - v)*dt + xi*sqrt(v_pos)*Z2*sqdt; 
            S=S + dS; v=fmax_generic(v + dv, asT(0.0,v)); 
            path.push_back(S);
        } 
        return path; 
    }
    
    template <typename T>
    /** @brief Simulate a Heston path from explicit variable values. */
    std::vector<T> simulate_path_Heston_from_vars(const OptionSpec& opt, const HestonSpec& H, const std::vector<double>& Z,
                                                const T& S0v, const T& rv, const T& kappav, const T& thetav, const T& xiv, const T& v0v, const T& rhov){
        const int n=opt.steps; 
        const double dt=opt.T/n; 
        const double sqdt=std::sqrt(dt);

        std::vector<T> path; 
        path.reserve(n); 
        T S=S0v; 
        T v=v0v; 
        for(int i=0;i<n;++i){ 
            double W1=Z[2*i], W2=Z[2*i+1]; 
            T Z1=asT(W1,S); 
            T Z2 = rhov*Z1 + sqrt( fmax_generic(asT(1.0,S) - rhov*rhov, asT(0.0,S)) ) * asT(W2,S);
            T v_pos=fmax_generic(v, asT(0.0,v)); 
            T dS=(rv - opt.q)*S*dt + sqrt(v_pos)*S*Z1*sqdt; 
            T dv=kappav*(thetav - v)*dt + xiv*sqrt(v_pos)*Z2*sqdt; 
            S=S + dS; v=fmax_generic(v + dv, asT(0.0,v)); 
            path.push_back(S);
        } 
        return path; 
    }
} 
#endif
