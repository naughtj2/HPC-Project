/**
 * @file aad.h
 * @brief Lightweight reverse-mode AD (tape-based) primitives.
 */
#ifndef PRICER_AAD_H
#define PRICER_AAD_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <type_traits>

namespace pricer {

    /**
     * @brief Computational tape holding nodes and adjoints for reverse-mode AD.
     */
    struct Tape { 
        /** @brief Operation kind encoded on the tape. */
        enum Op:uint8_t { CONST, VAR, ADD, SUB, MUL, DIV, EXP, LOG, SQRT, MAX }; 
        /** @brief Single tape node with operator and inputs.
         *  @param op Operator
         *  @param a  Left child index (or -1)
         *  @param b  Right child index (or -1)
         *  @param val Forward value cache
         */
        struct Node{Op op; int a,b; double val;}; 
        std::vector<Node> nodes; 
        std::vector<double> adj; 
        /** @brief Append a node to the tape and return its index. */
        int new_node(Op op,int a,int b,double val){ 
            nodes.push_back({op,a,b,val}); 
            return (int)nodes.size()-1; 
        } 
        /** @brief Erase all nodes and adjoints. */
        void clear(){ 
            nodes.clear(); 
            adj.clear(); 
        } 
    };

    /**
     * @brief Active double bound to a tape node.
     */
    struct ADouble { 
        int idx; 
        Tape* tape; 
        ADouble():idx(-1),tape(nullptr){} 
        ADouble(double v,Tape*t,bool is_var=false):tape(t){ 
            if(!t){idx=-1;return;} 
            idx=t->new_node(is_var?Tape::VAR:Tape::CONST,-1,-1,v);
        } 
        ADouble(int node_idx,Tape*t){
            idx=node_idx; 
            tape=t;
        }
        /** @brief Forward value of this node. */
        inline double val() const { 
            return tape->nodes[idx].val; 
        } 
    };

    /** @brief Create a constant on the given tape. */
    inline ADouble make_const(double v,Tape*t){
        return ADouble(v,t,false);
    } 
    /** @brief Create a variable (input) on the given tape. */
    inline ADouble make_var(double v,Tape*t){
        return ADouble(v,t,true);
    }
    /** @brief Addition of two active values. */
    ADouble operator+(const ADouble&x,const ADouble&y);
    /** @brief Subtraction of two active values. */
    ADouble operator-(const ADouble&x,const ADouble&y); 
    /** @brief Multiplication of two active values. */
    ADouble operator*(const ADouble&x,const ADouble&y); 
    /** @brief Division of two active values. */
    ADouble operator/(const ADouble&x,const ADouble&y); 
    inline /** @brief Subtraction of two active values. */
    ADouble operator-(const ADouble&x); 
    /** @brief Exponential function. */
    ADouble exp(const ADouble&x); 
    /** @brief Natural logarithm. */
    inline ADouble log(const ADouble&x); 
    /** @brief Square root. */
    ADouble sqrt(const ADouble&x); 
    /** @brief Pointwise maximum. */
    ADouble max(const ADouble&a,const ADouble&b); 
    
    /** @brief Reverse sweep to accumulate adjoints starting from node @p w. */
    void reverse_ad(Tape&t,int w);

    inline /** @brief Addition of two active values. */
    ADouble operator+(const ADouble&x,double c){
        return x+make_const(c,x.tape);
    } 
    inline /** @brief Addition of two active values. */
    ADouble operator+(double c,const ADouble&x){
        return make_const(c,x.tape)+x;
    } 
    inline /** @brief Subtraction of two active values. */
    ADouble operator-(const ADouble&x,double c){
        return x-make_const(c,x.tape);
    } 
    inline /** @brief Subtraction of two active values. */
    ADouble operator-(double c,const ADouble&x){
        return make_const(c,x.tape)-x;
    } 
    inline /** @brief Multiplication of two active values. */
    ADouble operator*(const ADouble&x,double c){
        return x*make_const(c,x.tape);
    } 
    inline /** @brief Multiplication of two active values. */
    ADouble operator*(double c,const ADouble&x){
        return make_const(c,x.tape)*x;
    } 
    inline /** @brief Division of two active values. */
    ADouble operator/(const ADouble&x,double c){
        return x/make_const(c,x.tape);
    } 
    inline /** @brief Division of two active values. */
    ADouble operator/(double c,const ADouble&x){
        return make_const(c,x.tape)/x;
    }

    /** @brief Overload-friendly fmax for doubles. */
    inline double fmax_generic(double a,double b){
        return std::max(a,b);
    } 
    /** @brief Overload-friendly fmax for active values. */
    inline ADouble fmax_generic(const ADouble&a,const ADouble&b){
        return max(a,b);
    }

    /** @brief Cast helper for templated math (generic type). */
    template <typename T> inline T asT(double v, const T&){ 
        return T(v); 
    } 
    /** @brief Cast helper: make an ADouble constant on the same tape as @p ref. */
    inline ADouble asT(double v, const ADouble& ref){ 
        return make_const(v, ref.tape); 
    }
} 
#endif
