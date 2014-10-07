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
  class CASADI_CORE_EXPORT SXFunction : public Function {

  public:
    /// Default constructor
    SXFunction();

    /// Expand an MXFunction
    explicit SXFunction(const MXFunction &f);

    /// Expand an Function
    explicit SXFunction(const Function &f);

    /// Multiple (matrix valued) input, multiple (matrix valued) output
    SXFunction(const std::vector< SX>& arg, const std::vector<SX>& res);

    /// Multiple (matrix valued) input, multiple (matrix valued) output
    SXFunction(const std::vector< SX>& arg, const IOSchemeVector< SX >& res);

    /// Multiple (matrix valued) input, multiple (matrix valued) output
    SXFunction(const IOSchemeVector< SX >& arg, const std::vector< SX>& res);

    /// Multiple (matrix valued) input, multiple (matrix valued) output
    SXFunction(const IOSchemeVector< SX >& arg, const IOSchemeVector< SX >& res);

#ifndef SWIG

    /// Multiple (vector valued) input, multiple (vector valued) output
    SXFunction(const std::vector< std::vector<SXElement> >& arg,
               const std::vector< std::vector<SXElement> >& res);

    /// Single (scalar/matrix/vector valued) input, single (scalar/matrix/vector valued) output
    SXFunction(const SX& arg, const SX& res);

    /// Multiple (vector valued) input, single (scalar/vector/matrix valued) output
    SXFunction(const std::vector< std::vector<SXElement> >& arg, const SX& res);

    /// Multiple (matrix valued) input, single (scalar/vector/matrix valued) output
    SXFunction(const std::vector< SX>& arg, const SX& res);

    /// Single (scalar/vector/matrix valued) input, multiple (vector valued) output
    SXFunction(const SX& arg, const std::vector< std::vector<SXElement> >& res);

    /// Single (scalar/vector/matrix valued) input, multiple (matrix valued) output
    SXFunction(const SX& arg, const std::vector< SX>& res);
#endif // SWIG

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
    { return jac(inputSchemeEntry(iname), oind, compact, symmetric); }
    SX jac(int iind, const std::string& oname, bool compact=false, bool symmetric=false)
    { return jac(iind, outputSchemeEntry(oname), compact, symmetric); }
    SX jac(const std::string& iname, const std::string& oname,
           bool compact=false, bool symmetric=false)
    { return jac(inputSchemeEntry(iname), outputSchemeEntry(oname), compact, symmetric); }
    ///@}

    ///@{
    /// Gradient via source code transformation
    SX grad(int iind=0, int oind=0);
    SX grad(const std::string& iname, int oind=0) { return grad(inputSchemeEntry(iname), oind); }
    SX grad(int iind, const std::string& oname) { return grad(iind, outputSchemeEntry(oname)); }
    SX grad(const std::string& iname, const std::string& oname)
    { return grad(inputSchemeEntry(iname), outputSchemeEntry(oname)); }
    ///@}

    ///@{
    /// Tangent via source code transformation
    SX tang(int iind=0, int oind=0);
    SX tang(const std::string& iname, int oind=0) { return tang(inputSchemeEntry(iname), oind); }
    SX tang(int iind, const std::string& oname) { return tang(iind, outputSchemeEntry(oname)); }
    SX tang(const std::string& iname, const std::string& oname)
    { return tang(inputSchemeEntry(iname), outputSchemeEntry(oname)); }
    ///@}

    ///@{
    /// Hessian (forward over adjoint) via source code transformation
    SX hess(int iind=0, int oind=0);
    SX hess(const std::string& iname, int oind=0) { return hess(inputSchemeEntry(iname), oind); }
    SX hess(int iind, const std::string& oname) { return hess(iind, outputSchemeEntry(oname)); }
    SX hess(const std::string& iname, const std::string& oname)
    { return hess(inputSchemeEntry(iname), outputSchemeEntry(oname)); }
    ///@}

    /** \brief Get function input */
    const SX& inputExpr(int iind) const;
    const SX& inputExpr(const std::string& iname) const
    { return inputExpr(inputSchemeEntry(iname)); }

    /** \brief Get function output */
    const SX& outputExpr(int oind) const;
    const SX& outputExpr(const std::string& oname) const
    { return outputExpr(outputSchemeEntry(oname)); }

    /** \brief Get all function inputs */
    const std::vector<SX>& inputExpr() const;

    /** \brief Get all function outputs */
    const std::vector<SX> & outputExpr() const;

/// \cond INTERNAL
#ifndef SWIG
    /** \brief Access the algorithm directly */
    const std::vector<ScalarAtomic>& algorithm() const;
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
    std::vector<SXElement> getFree() const;

    /** \brief Get the corresponding matrix type */
    typedef SX MatType;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
  };

} // namespace casadi

#endif // CASADI_SX_FUNCTION_HPP
