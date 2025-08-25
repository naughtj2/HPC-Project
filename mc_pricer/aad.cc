/**
 * @file aad.cc
 * @brief Operators and reverse pass for the lightweight reverse-mode AD types.
 */
#include "aad.h"

namespace pricer {

    inline ADouble _mk(Tape::Op op, const ADouble& x, const ADouble& y, double v){
        int id = x.tape->new_node(op, x.idx, y.idx, v);
        return ADouble(id, x.tape);
    }

    /**
 * @brief Overloaded arithmetic operators create new nodes on the tape and cache forward values.
 */
    ADouble operator+(const ADouble& x, const ADouble& y) { return _mk(Tape::ADD, x, y, x.val()+y.val()); }
    ADouble operator-(const ADouble& x, const ADouble& y) { return _mk(Tape::SUB, x, y, x.val()-y.val()); }
    ADouble operator*(const ADouble& x, const ADouble& y) { return _mk(Tape::MUL, x, y, x.val()*y.val()); }
    ADouble operator/(const ADouble& x, const ADouble& y) { return _mk(Tape::DIV, x, y, x.val()/y.val()); }

    ADouble operator-(const ADouble& x) {
        Tape* t = x.tape; ADouble zero = make_const(0.0, t);
        int id = t->new_node(Tape::SUB, zero.idx, x.idx, -x.val());
        return ADouble(id, t);
    }

    /** @brief Exponential of an active value. */
    ADouble exp(const ADouble& x) {
        double v = std::exp(x.val());
        int id = x.tape->new_node(Tape::EXP, x.idx, -1, v);
        return ADouble(id, x.tape);
    }
    /** @brief Natural logarithm of an active value. */
    ADouble log(const ADouble& x) {
        double v = std::log(x.val());
        int id = x.tape->new_node(Tape::LOG, x.idx, -1, v);
        return ADouble(id, x.tape);
    }
    /** @brief Square root of an active value. */
    ADouble sqrt(const ADouble& x) {
        double v = std::sqrt(x.val());
        int id = x.tape->new_node(Tape::SQRT, x.idx, -1, v);
        return ADouble(id, x.tape);
    }
    /** @brief Pointwise maximum of two active values. */
    ADouble max(const ADouble& x, const ADouble& y) {
        double v = std::max(x.val(), y.val());
        int id = x.tape->new_node(Tape::MAX, x.idx, y.idx, v);
        return ADouble(id, x.tape);
    }

    /**
 * @brief Reverse-mode sweep accumulating adjoints on @p tape starting from output node @p w.
 * @param tape The computational tape
 * @param w Index of the scalar output to differentiate
 */
    void reverse_ad(Tape& tape, int w) {
        const size_t N = tape.nodes.size();
        tape.adj.assign(N, 0.0);
        tape.adj[w] = 1.0;
        for (int i = (int)N - 1; i >= 0; --i) {
            double bar = tape.adj[i]; if (bar == 0.0) continue;
            const auto& n = tape.nodes[i];
            switch (n.op) {
                case Tape::ADD:
                    tape.adj[n.a] += bar;
                    tape.adj[n.b] += bar;
                    break;
                case Tape::SUB:
                    tape.adj[n.a] += bar;
                    tape.adj[n.b] -= bar;
                    break;
                case Tape::MUL: {
                    double aval = tape.nodes[n.a].val;
                    double bval = tape.nodes[n.b].val;
                    tape.adj[n.a] += bar * bval;
                    tape.adj[n.b] += bar * aval;
                    break; }
                case Tape::DIV: {
                    double aval = tape.nodes[n.a].val;
                    double bval = tape.nodes[n.b].val;
                    tape.adj[n.a] += bar / bval;
                    tape.adj[n.b] -= bar * aval / (bval * bval);
                    break; }
                case Tape::EXP: {
                    tape.adj[n.a] += bar * n.val;
                    break; }
                case Tape::LOG: {
                    double aval = tape.nodes[n.a].val;
                    tape.adj[n.a] += bar / aval;
                    break; }
                case Tape::SQRT: {
                    if (n.val > 0.0)
                        tape.adj[n.a] += bar * (0.5 / n.val);
                    break; }
                case Tape::MAX: {
                    double aval = tape.nodes[n.a].val;
                    double bval = tape.nodes[n.b].val;
                    if (aval > bval) tape.adj[n.a] += bar; else if (bval > aval) tape.adj[n.b] += bar;
                    break; }
                case Tape::CONST:
                case Tape::VAR:
                default: break;
            }
        }
    }

} 
