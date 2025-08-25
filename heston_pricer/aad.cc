/**
 * @file aad.cc
 * @brief Arithmetic operators and reverse sweep for AAD types.
 */
#include "aad.h"

namespace pricer {

    /**
     * @brief Internal helper to allocate a node and wrap it as @ref ADouble.
     */
    inline ADouble _mk(Tape::Op op, const ADouble& x, const ADouble& y, double v){ 
        int id = x.tape->new_node(op, x.idx, y.idx, v); 
        return ADouble(id, x.tape); 
        }

    /**
     * @brief Overloaded arithmetic creates new tape nodes and caches forward values.
     */
    ADouble operator+(const ADouble& x, const ADouble& y){ 
        return _mk(Tape::ADD, x, y, x.val()+y.val()); 
    }
    ADouble operator-(const ADouble& x, const ADouble& y){ 
        return _mk(Tape::SUB, x, y, x.val()-y.val()); 
    }
    ADouble operator*(const ADouble& x, const ADouble& y){ 
        return _mk(Tape::MUL, x, y, x.val()*y.val()); 
    }
    ADouble operator/(const ADouble& x, const ADouble& y){ 
        return _mk(Tape::DIV, x, y, x.val()/y.val());
     }
    ADouble operator-(const ADouble& x){ 
        Tape* t=x.tape; 
        ADouble zero=make_const(0.0,t); 
        int id=t->new_node(Tape::SUB, zero.idx, x.idx, -x.val()); 
        return ADouble(id,t);
    }
    ADouble exp(const ADouble& x){ 
        double v=std::exp(x.val()); 
        int id=x.tape->new_node(Tape::EXP, x.idx, -1, v); 
        return ADouble(id, x.tape); 
    }
    ADouble log(const ADouble& x){ 
        double v=std::log(x.val()); 
        int id=x.tape->new_node(Tape::LOG, x.idx, -1, v); 
        return ADouble(id, x.tape); 
    }
    ADouble sqrt(const ADouble& x){ double v=std::sqrt(x.val()); 
        int id=x.tape->new_node(Tape::SQRT, x.idx, -1, v); 
        return ADouble(id, x.tape); 
    }
    ADouble max(const ADouble& a, const ADouble& b){ 
        double v=std::max(a.val(), b.val()); 
        int id=a.tape->new_node(Tape::MAX, a.idx, b.idx, v); 
        return ADouble(id, a.tape); 
    }
    /**
    * @brief Reverse-mode accumulation of adjoints on @p t starting from output @p w.
    */
    void reverse_ad(Tape&t,int w){ 
        size_t N=t.nodes.size(); 
        t.adj.assign(N,0.0); 
        t.adj[w]=1.0; 
        for(int i=(int)N-1;i>=0;--i){ 
            double bar=t.adj[i]; 
            if(bar==0.0) continue; 
            const auto& n=t.nodes[i]; 
            switch(n.op){ 
                case Tape::ADD: 
                t.adj[n.a]+=bar; 
                t.adj[n.b]+=bar; 
                break; 
                case Tape::SUB: t.adj[n.a]+=bar; 
                t.adj[n.b]-=bar; 
                break; 
                case Tape::MUL:{ 
                    double av=t.nodes[n.a].val, bv=t.nodes[n.b].val; 
                    t.adj[n.a]+=bar*bv; 
                    t.adj[n.b]+=bar*av; 
                    break;
                } 
                case Tape::DIV:{ 
                    double av=t.nodes[n.a].val, bv=t.nodes[n.b].val; 
                    t.adj[n.a]+=bar/bv; 
                    t.adj[n.b]-=bar*av/(bv*bv); 
                    break;
                } 
                case Tape::EXP: t.adj[n.a]+=bar*n.val; 
                break; 
                case Tape::LOG: t.adj[n.a]+=bar/t.nodes[n.a].val; 
                break; 
                case Tape::SQRT: if(n.val>0.0) t.adj[n.a]+=bar*(0.5/n.val); 
                break; 
                case Tape::MAX:{ 
                    double av=t.nodes[n.a].val, bv=t.nodes[n.b].val; 
                    if(av>bv) t.adj[n.a]+=bar; 
                    else if(bv>av) t.adj[n.b]+=bar; 
                    break;
                } 
                default: break; 
            } 
        } 
    }
} 