/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_IMPLICIT_FUNCTION_INTERNAL_HPP
#define CASADI_IMPLICIT_FUNCTION_INTERNAL_HPP

#include "implicit_function.hpp"
#include "function_internal.hpp"
#include "linear_solver.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

  /// Internal class
  class CASADI_EXPORT
  ImplicitFunctionInternal : public FunctionInternal,
                             public PluginInterface<ImplicitFunctionInternal> {
  public:
    /** \brief Constructor
     *
     * \param f   Function mapping from (n+1) inputs to 1 output.
     */
    ImplicitFunctionInternal(const Function& f);

    /// Destructor
    virtual ~ImplicitFunctionInternal() = 0;

    /** \brief  Clone */
    virtual ImplicitFunctionInternal* clone() const=0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new ImplicitFunctionInternal */
    virtual ImplicitFunctionInternal* create(const Function& f) const = 0;

    /// Initialize
    virtual void init();

    /** \brief  Propagate the sparsity pattern through a set of
     * directional derivatives forward or backward */
    virtual void spEvaluate(bool fwd);

    /// Is the class able to propagate seeds through the algorithm?
    virtual bool spCanEvaluate(bool fwd) { return true;}

    /// Solve the system of equations and calculate derivatives
    virtual void evaluate();

    /// Solve the nonlinear system of equations
    virtual void solveNonLinear() = 0;

    ///@{
    /** \brief Generate a function that calculates nfwd forward derivatives */
    virtual Function getDerivativeFwd(int nfwd);
    virtual bool hasDerivativeFwd() const { return true;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates nadj adjoint derivatives */
    virtual Function getDerivativeAdj(int nadj);
    virtual bool hasDerivativeAdj() const { return true;}
    ///@}

    /** \brief Create call to (cached) derivative function, forward mode  */
    virtual void callFwd(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens,
                         bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    virtual void callAdj(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens,
                         bool always_inline, bool never_inline);

    /// Number of equations
    int n_;

    /// The function f(z, x1, x2, ..., xn) == 0
    Function f_;

    /// Jacobian of f with respect to z
    Function jac_;

    /// Linear solver
    LinearSolver linsol_;

    /// Factorization up-to-date?
    bool fact_up_to_date_;

    /// Constraints on decision variables
    std::vector<int> u_c_;

    /// Indices of the input and output that correspond to the actual root-finding
    int iin_, iout_;

    // Creator function for internal class
    typedef ImplicitFunctionInternal* (*Creator)(const Function& f);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Short name
    static std::string shortname() { return "nlsol";}

    /// Infix
    static const std::string infix_;
  };



} // namespace casadi
/// \endcond

#endif // CASADI_IMPLICIT_FUNCTION_INTERNAL_HPP

