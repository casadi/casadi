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


#ifndef CASADI_SX_FUNCTION_HPP
#define CASADI_SX_FUNCTION_HPP

#include "function.hpp"
#ifdef USE_CXX11
#include <initializer_list>
#endif // USE_CXX11

namespace casadi {

/// \cond INTERNAL
#ifndef SWIG

  /** \brief  An atomic operation for the SXElement virtual machine */
  struct ScalarAtomic {
    int op;     /// Operator index
    int i0;
    union {
      double d;
      struct { int i1, i2; };
    };
  };

#endif // SWIG

/// \endcond

  /// Forward declaration of internal class
  class SXFunctionInternal;

  /// Forward declaration of MXFunction
  class MXFunction;

  /**   \brief Dynamically created function that can be expanded into a series of scalar operations.

        \author Joel Andersson
        \date 2010-2013
  */
  class CASADI_EXPORT SXFunction : public Function {

  public:
    /// Default constructor
    SXFunction();

    /// Expand an MXFunction
    explicit SXFunction(const MXFunction &f);

    /// Expand an Function
    explicit SXFunction(const Function &f);

    /** \brief Construct from vectors (new syntax, includes initialization) */
    SXFunction(const std::string& name, const std::vector<SX>& arg,
               const std::vector<SX>& res, const Dict& opts=Dict());

    /** \brief Construct from vectors (new syntax, includes initialization) */
    SXFunction(const std::string& name, const std::pair< SXDict, std::vector<std::string> >& arg,
               const std::vector<SX>& res, const Dict& opts=Dict());

    /** \brief Construct from vectors (new syntax, includes initialization) */
    SXFunction(const std::string& name, const std::vector<SX>& arg,
               const std::pair< SXDict, std::vector<std::string> >& res, const Dict& opts=Dict());

    /** \brief Construct from vectors (new syntax, includes initialization) */
    SXFunction(const std::string& name, const std::pair< SXDict, std::vector<std::string> >& arg,
               const std::pair< SXDict, std::vector<std::string> >& res, const Dict& opts=Dict());
#ifndef SWIG
#ifdef USE_CXX11
    /** \brief Construct from initializer lists (new syntax, includes initialization) */
    SXFunction(const std::string& name,
               std::initializer_list<SX> arg,
               std::initializer_list<SX> res,
               const Dict& opts=Dict());

    /** \brief Construct from vector & nitializer list (new syntax, includes initialization) */
    SXFunction(const std::string& name,
               std::vector<SX> arg,
               std::initializer_list<SX> res,
               const Dict& opts=Dict());

    /** \brief Construct from initializer list & vector (new syntax, includes initialization) */
    SXFunction(const std::string& name,
               std::initializer_list<SX> arg,
               std::vector<SX> res,
               const Dict& opts=Dict());
#endif // USE_CXX11
#endif // SWIG

#ifdef WITH_DEPRECATED_FEATURES
    /// [DEPRECATED] Multiple input, multiple output, no initialization
    SXFunction(const std::vector<SX>& arg,
               const std::vector<SX>& res);

    /// [DEPRECATED] Multiple input, multiple  output, no initialization
    SXFunction(const std::vector<SX>& arg,
               const std::pair< SXDict, std::vector<std::string> >& res);

    /// [DEPRECATED] Multiple input, multiple output, no initialization
    SXFunction(const std::pair< SXDict, std::vector<std::string> >& arg,
               const std::vector<SX>& res);

    /// [DEPRECATED] Multiple input, multiple output, no initialization
    SXFunction(const std::pair< SXDict, std::vector<std::string> >& arg,
               const std::pair< SXDict, std::vector<std::string> >& res);
#endif // WITH_DEPRECATED_FEATURES

/// \cond INTERNAL
    /// Access functions of the node
    SXFunctionInternal* operator->();

    /// Const access functions of the node
    const SXFunctionInternal* operator->() const;
/// \endcond

    ///@{
    /** \brief Jacobian via source code transformation
     *
     * \see casadi::Jacobian for an AD approach
     */
    SX jac(int iind=0, int oind=0, bool compact=false, bool symmetric=false);
    SX jac(const std::string& iname, int oind=0, bool compact=false, bool symmetric=false)
    { return jac(inputIndex(iname), oind, compact, symmetric); }
    SX jac(int iind, const std::string& oname, bool compact=false, bool symmetric=false)
    { return jac(iind, outputIndex(oname), compact, symmetric); }
    SX jac(const std::string& iname, const std::string& oname,
           bool compact=false, bool symmetric=false)
    { return jac(inputIndex(iname), outputIndex(oname), compact, symmetric); }
    ///@}

    ///@{
    /// Gradient via source code transformation
    SX grad(int iind=0, int oind=0);
    SX grad(const std::string& iname, int oind=0) { return grad(inputIndex(iname), oind); }
    SX grad(int iind, const std::string& oname) { return grad(iind, outputIndex(oname)); }
    SX grad(const std::string& iname, const std::string& oname)
    { return grad(inputIndex(iname), outputIndex(oname)); }
    ///@}

    ///@{
    /// Tangent via source code transformation
    SX tang(int iind=0, int oind=0);
    SX tang(const std::string& iname, int oind=0) { return tang(inputIndex(iname), oind); }
    SX tang(int iind, const std::string& oname) { return tang(iind, outputIndex(oname)); }
    SX tang(const std::string& iname, const std::string& oname)
    { return tang(inputIndex(iname), outputIndex(oname)); }
    ///@}

    ///@{
    /// Hessian (forward over adjoint) via source code transformation
    SX hess(int iind=0, int oind=0);
    SX hess(const std::string& iname, int oind=0) { return hess(inputIndex(iname), oind); }
    SX hess(int iind, const std::string& oname) { return hess(iind, outputIndex(oname)); }
    SX hess(const std::string& iname, const std::string& oname) {
      return hess(inputIndex(iname), outputIndex(oname));
    }
    ///@}

    /** \brief Get function input */
    const SX inputExpr(int iind) const;
    const SX inputExpr(const std::string& iname) const {
      return inputExpr(inputIndex(iname));
    }

    /** \brief Get function output */
    const SX outputExpr(int oind) const;
    const SX outputExpr(const std::string& oname) const {
      return outputExpr(outputIndex(oname));
    }

    /** \brief Get all function inputs */
    const std::vector<SX> inputExpr() const;

    /** \brief Get all function outputs */
    const std::vector<SX> outputExpr() const;

/// \cond INTERNAL
#ifndef SWIG
    /** \brief Access the algorithm directly */
    const std::vector<ScalarAtomic>& algorithm() const;

    /** \brief Called from constructor */
    void construct(const std::string& name, const std::vector<SX>& arg,
                   const std::vector<SX>& res, const Dict& opts,
                   const std::vector<std::string>& ischeme=std::vector<std::string>(),
                   const std::vector<std::string>& oscheme=std::vector<std::string>());
#endif // SWIG
/// \endcond

    /** \brief Get the number of atomic operations */
    int getAlgorithmSize() const { return algorithm().size();}

    /** \brief Get the length of the work vector */
    int getWorkSize() const;

    /** \brief Get an atomic operation operator index */
    int getAtomicOperation(int k) const { return algorithm().at(k).op;}

    /** \brief Get the (integer) input arguments of an atomic operation */
    std::pair<int, int> getAtomicInput(int k) const {
      const ScalarAtomic& atomic = algorithm().at(k);
      return std::pair<int, int>(atomic.i1, atomic.i2);}

    /** \brief Get the floating point output argument of an atomic operation */
    double getAtomicInputReal(int k) const { return algorithm().at(k).d;}

    /** \brief Get the (integer) output argument of an atomic operation */
    int getAtomicOutput(int k) const { return algorithm().at(k).i0;}

    /** \brief Number of nodes in the algorithm */
    int countNodes() const;

    /** \brief Clear the function from its symbolic representation, to free up memory,
     * no symbolic evaluations are possible after this */
    void clearSymbolic();

    /** \brief Get all the free variables of the function */
    SX getFree() const;

    /** \brief Get the corresponding matrix type */
    typedef SX MatType;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
  };

} // namespace casadi

#endif // CASADI_SX_FUNCTION_HPP
