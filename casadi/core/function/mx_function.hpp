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


#ifndef CASADI_MX_FUNCTION_HPP
#define CASADI_MX_FUNCTION_HPP

#include <set>
#include <iostream>

#include "../mx/mx.hpp"
#include "sx_function.hpp"

/// \cond INTERNAL

namespace casadi {

/** \brief  An element of the algorithm, namely an MX node */
struct MXAlgEl {
  /// Operator index
  int op;

  /// Data associated with the operation
  MX data;

  /// Work vector indices of the arguments
  std::vector<int> arg;

  /// Work vector indices of the results
  std::vector<int> res;
};

} // namespace casadi

#ifdef SWIG
// Template instantiation
%template(MXAlgElVector) std::vector<casadi::MXAlgEl>;
#endif // SWIG

/// \endcond

namespace casadi {

/** \brief  Forward declaration of internal class */
class MXFunctionInternal;

  /** \brief  General function mapping from/to MX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT MXFunction : public Function {
public:

  /** \brief  Default constructor */
  MXFunction();

  /** \brief  Attempt to form an MXFunction out of an Function */
  explicit MXFunction(const Function& function);

#ifndef SWIG
  /** \brief  Single input, single output */
  MXFunction(const MX& input, const MX& output);

  /** \brief  Single input, multiple output */
  MXFunction(const MX& input, const std::vector<MX>& output);

  /** \brief  Multiple input, single output */
  MXFunction(const std::vector<MX>& input, const MX& output);
#endif // SWIG

  /** \brief  Multiple input, multiple output*/
  MXFunction(const std::vector<MX>& input, const std::vector<MX>& output);

  /** \brief  Multiple input, multiple output*/
  MXFunction(const std::vector<MX>& input, const IOSchemeVector< MX >& output);

  /** \brief  Multiple input, multiple output*/
  MXFunction(const IOSchemeVector< MX >& input, const std::vector<MX>& output);

  /** \brief  Multiple input, multiple output*/
  MXFunction(const IOSchemeVector< MX >& input, const IOSchemeVector< MX >& output);

  /// \cond INTERNAL
  /** \brief  Access functions of the node */
  MXFunctionInternal* operator->();

  /** \brief  Const access functions of the node */
  const MXFunctionInternal* operator->() const;
  /// \endcond

  /** \brief Get function input */
  const MX& inputExpr(int ind) const;
  const MX& inputExpr(const std::string & iname) const
  { return inputExpr(inputSchemeEntry(iname)); }

  /** \brief Get function output */
  const MX& outputExpr(int ind) const;
  const MX& outputExpr(const std::string & oname) const
  { return outputExpr(outputSchemeEntry(oname)); }

  /** \brief Get all function inputs */
  const std::vector<MX>& inputExpr() const;

  /** \brief Get all function outputs */
  const std::vector<MX> & outputExpr() const;

/// \cond INTERNAL
#ifndef SWIG
  /** \brief Access the algorithm directly */
  const std::vector<MXAlgEl>& algorithm() const;
#endif // SWIG
/// \endcond

  /** \brief Get the number of atomic operations */
  int getAlgorithmSize() const { return algorithm().size();}

  /** \brief Get the length of the work vector */
  int getWorkSize() const;

  /** \brief Number of nodes in the algorithm */
  int countNodes() const;

  /// Check if a particular cast is allowed
  static bool testCast(const SharedObjectNode* ptr);

  ///@{
  /** \brief Jacobian via source code transformation */
  MX jac(int iind=0, int oind=0, bool compact=false, bool symmetric=false);
  MX jac(const std::string & iname, int oind=0, bool compact=false, bool symmetric=false)
  { return jac(inputSchemeEntry(iname), oind, compact, symmetric); }
  MX jac(int iind, const std::string & oname, bool compact=false, bool symmetric=false)
  { return jac(iind, outputSchemeEntry(oname), compact, symmetric); }
  MX jac(const std::string & iname, const std::string & oname,
         bool compact=false, bool symmetric=false)
  { return jac(inputSchemeEntry(iname), outputSchemeEntry(oname), compact, symmetric); }
  ///@}

  ///@{
  /** \brief Gradient via source code transformation */
  MX grad(int iind=0, int oind=0);
  MX grad(const std::string & iname, int oind=0) { return grad(inputSchemeEntry(iname), oind); }
  MX grad(int iind, const std::string & oname) { return grad(iind, outputSchemeEntry(oname)); }
  MX grad(const std::string & iname, const std::string & oname)
  { return grad(inputSchemeEntry(iname), outputSchemeEntry(oname)); }
  ///@}

  ///@{
  /** \brief Tangent via source code transformation */
  MX tang(int iind=0, int oind=0);
  MX tang(const std::string & iname, int oind=0) { return tang(inputSchemeEntry(iname), oind); }
  MX tang(int iind, const std::string & oname) { return tang(iind, outputSchemeEntry(oname)); }
  MX tang(const std::string & iname, const std::string & oname)
  { return tang(inputSchemeEntry(iname), outputSchemeEntry(oname)); }
  ///@}

  /** \brief Expand the matrix valued graph into a scalar valued graph */
  SXFunction expand(const std::vector<SX>& inputv = std::vector<SX>());

  /** \brief Get all the free variables of the function */
  std::vector<MX> getFree() const;

/// \cond INTERNAL
#ifndef SWIG
  /** \brief Extract the functions needed for the Lifted Newton method */
  void generateLiftingFunctions(MXFunction& vdef_fcn, MXFunction& vinit_fcn);
#else
  /** \brief Extract the functions needed for the Lifted Newton method */
  void generateLiftingFunctions(MXFunction& OUTPUT, MXFunction& OUTPUT);
#endif
/// \endcond

  /** \brief Get the corresponding matrix type */
  typedef MX MatType;
};

} // namespace casadi


#endif // CASADI_MX_FUNCTION_HPP

