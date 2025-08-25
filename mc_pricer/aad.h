/**
 * @file aad.h
 * @brief Lightweight reverse-mode AD (tape-based) types and operators.
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
 * @brief Computational tape storing operation nodes and adjoints for reverse-mode AD.
 */
    struct Tape {
        /** @brief Operation code for a node on the tape. */
        enum Op : uint8_t { CONST, VAR, ADD, SUB, MUL, DIV, EXP, LOG, SQRT, MAX };
        /** @brief A single node on the tape.
         *  @param op Operation type
         *  @param a Index of left child (or -1)
         *  @param b Index of right child (or -1)
         *  @param val Cached forward value
         */
        struct Node { Op op; int a, b; double val; };
        std::vector<Node> nodes;
        std::vector<double> adj;
        /** @brief Append a node to the tape.
         *  @return Index of the newly created node.
         */
        inline int new_node(Op op, int a, int b, double val) {
            nodes.push_back({op, a, b, val});
            return (int)nodes.size() - 1;
        }
        /** @brief Reset the tape contents. */
        void clear() { nodes.clear(); adj.clear(); }
    };

    /**
 * @brief Active double value tied to a tape node.
 *
 * Wraps an index into a @ref Tape and exposes the cached forward value.
 */
    struct ADouble {
        int idx;
        Tape* tape;
        ADouble() : idx(-1), tape(nullptr) {}
        ADouble(double v, Tape* t, bool is_var=false) : tape(t) {
            if (!t) { idx = -1; return; }
            idx = t->new_node(is_var ? Tape::VAR : Tape::CONST, -1, -1, v);
        }
        ADouble(int node_idx, Tape* t) { idx = node_idx; tape = t; }
        /** @brief Forward value cached at this node. */
        inline double val() const { return tape->nodes[idx].val; }
    };

    /** @brief Create a constant node on the tape. */
    inline ADouble make_const(double v, Tape* t) { return ADouble(v, t, false); }
    /** @brief Create a variable (input) node on the tape. */
    inline ADouble make_var(double v, Tape* t)   { return ADouble(v, t, true); }

    /** @brief Add two active values. */
    ADouble operator+(const ADouble& x, const ADouble& y);
    /** @brief Subtract two active values. */
    ADouble operator-(const ADouble& x, const ADouble& y);
    /** @brief Multiply two active values. */
    ADouble operator*(const ADouble& x, const ADouble& y);
    /** @brief Divide two active values. */
    ADouble operator/(const ADouble& x, const ADouble& y);
    inline /** @brief Subtract two active values. */
    ADouble operator-(const ADouble& x);
    /** @brief Exponential function. */
    ADouble exp(const ADouble& x);
    /** @brief Natural logarithm. */
    inline ADouble log(const ADouble& x);
    /** @brief Square root. */
    ADouble sqrt(const ADouble& x);
    /** @brief Pointwise maximum. */
    ADouble max(const ADouble& x, const ADouble& y);
    /** @brief Perform reverse-mode sweep to fill tape adjoints starting from output index @p w. */
    void reverse_ad(Tape& tape, int w);

    inline /** @brief Add two active values. */
    ADouble operator+(const ADouble& x, double c) { return x + make_const(c, x.tape); }
    inline /** @brief Add two active values. */
    ADouble operator+(double c, const ADouble& x) { return make_const(c, x.tape) + x; }
    inline /** @brief Subtract two active values. */
    ADouble operator-(const ADouble& x, double c) { return x - make_const(c, x.tape); }
    inline /** @brief Subtract two active values. */
    ADouble operator-(double c, const ADouble& x) { return make_const(c, x.tape) - x; }
    inline /** @brief Multiply two active values. */
    ADouble operator*(const ADouble& x, double c) { return x * make_const(c, x.tape); }
    inline /** @brief Multiply two active values. */
    ADouble operator*(double c, const ADouble& x) { return make_const(c, x.tape) * x; }
    inline /** @brief Divide two active values. */
    ADouble operator/(const ADouble& x, double c) { return x / make_const(c, x.tape); }
    inline /** @brief Divide two active values. */
    ADouble operator/(double c, const ADouble& x) { return make_const(c, x.tape) / x; }

    /** @brief Overload-friendly fmax for double. */
    inline double fmax_generic(double a, double b) { return std::max(a,b); }
    /** @brief Overload-friendly fmax for active values. */
    inline ADouble fmax_generic(const ADouble& a, const ADouble& b) { return max(a,b); }

    template <typename T>
    /** @brief Convert double to type @p T (generic). */
    inline T asT(double v, const T&) { return T(v); }
    /** @brief Convert double into an @ref ADouble constant on the same tape as @p ref. */
    inline ADouble asT(double v, const ADouble& ref) { return make_const(v, ref.tape); }

} 

#endif
