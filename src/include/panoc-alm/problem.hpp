#pragma once

#include "box.hpp"

#include <functional>

namespace pa {

/**
 * @class Problem
 * @brief   Problem description for minimization problems.
 * 
 * @f[ \begin{aligned}
 *  & \underset{x}{\text{minimize}}
 *  & & f(x) &&&& f : \mathbb{R}^n \rightarrow \mathbb{R} \\
 *  & \text{subject to}
 *  & & \underline{x} \le x \le \overline{x} \\
 *  &&& \underline{z} \le g(x) \le \overline{z} &&&& g :
 *  \mathbb{R}^n \rightarrow \mathbb{R}^m
 * \end{aligned} @f]
 */
struct Problem {
    unsigned int n; ///< Number of decision variables, dimension of x
    unsigned int m; ///< Number of constraints, dimension of g(x) and z
    Box C;          ///< Constraints of the decision variables, @f$ x \in C @f$
    Box D;          ///< Other constraints, @f$ g(x) \in D @f$

    /// Signature of the function that evaluates the cost @f$ f(x) @f$
    /// @param  x
    ///   [in]  Decision variable: @f$ x @f$
    using f_sig = real_t(const vec &x);
    /// Signature of the function that evaluates the gradient of the cost
    /// function @f$ \nabla f(x) @f$
    /// @param  x
    ///   [in]  Decision variable: @f$ x @f$
    /// @param  grad_fx
    ///  [out]  Gradient of cost function: @f$ \nabla f(x) @f$
    using grad_f_sig = void(const vec &x, vec &grad_fx);
    /// Signature of the function that evaluates the constraints @f$ g(x) @f$
    /// @param  x
    ///   [in]  Decision variable: @f$ x @f$
    /// @param  gx
    ///  [out]  Value of the constraints: @f$ g(x) @f$
    using g_sig = void(const vec &x, vec &gx);
    /// Signature of the function that evaluates the gradient of the constraints
    /// times a vector
    /// @f$ \nabla g(x)\ y @f$
    /// @param  x
    ///   [in]  Decision variable: @f$ x @f$
    /// @param  y
    ///   [in]  Vector to multiply the gradient by: @f$ y @f$
    /// @param  grad_gxy
    ///  [out]  Gradient of the constraints: @f$ \nabla g(x)\ y @f$
    using grad_g_sig = void(const vec &x, const vec &y, vec &grad_gxy);

    /// Cost function @f$ f(x) @f$
    std::function<f_sig> f;
    /// Gradient of the cost function @f$ \nabla f(x) @f$
    std::function<grad_f_sig> grad_f;
    /// Constraint function @f$ g(x) @f$
    std::function<g_sig> g;
    /// Gradient of the constraint function times vector @f$ \nabla g(x)\ y @f$
    std::function<grad_g_sig> grad_g;
};

} // namespace pa
