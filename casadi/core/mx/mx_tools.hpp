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


#ifndef CASADI_MX_TOOLS_HPP
#define CASADI_MX_TOOLS_HPP

#include "mx.hpp"

#include "../matrix/generic_expression_tools.hpp"
#include "../function/linear_solver.hpp"

namespace casadi {

/**
\ingroup expression_tools
@{
*/
#if !defined(SWIG) || !defined(SWIGMATLAB)

  /** \brief Branching on MX nodes
      Ternary operator, "cond ? if_true : if_false"
  */
  inline MX if_else(const MX &cond, const MX &if_true, const MX &if_false) {
    return cond.zz_if_else(if_true, if_false);
  }

  /** \brief  Simplify an expression */
  inline void simplify(MX& ex) { ex.zz_simplify();}

  /** \brief  Create a parent MX on which all given MX's will depend.

      In some sense, this function is the inverse of

      \param deps  Must all be symbolic matrices.
  */
  inline MX createParent(std::vector<MX> &deps) { return MX::zz_createParent(deps);}

  /** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
   */
  inline MX createParent(const std::vector<MX> &deps,
                         std::vector<MX>& SWIG_OUTPUT(children)) {
    return MX::zz_createParent(deps, children);
  }

  /** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
   */
  inline MX createParent(const std::vector<Sparsity> &deps,
                         std::vector<MX>& SWIG_OUTPUT(children)) {
    return MX::zz_createParent(deps, children);
  }

  /** Count number of nodes */
  inline int countNodes(const MX& A) { return A.zz_countNodes();}

  /** \brief Get a string representation for a binary MX, using custom arguments */
  inline std::string
  getOperatorRepresentation(const MX& xb, const std::vector<std::string>& args) {
    return xb.zz_getOperatorRepresentation(args);
  }

  /** \brief  Substitute variable v with expression vdef in an expression ex */
  inline MX substitute(const MX &ex, const MX& v, const MX& vdef) {
    return ex.zz_substitute(v, vdef);
  }

  /** \brief  Substitute variable var with expression expr in multiple expressions */
  inline std::vector<MX>
  substitute(const std::vector<MX> &ex, const std::vector<MX> &v, const std::vector<MX> &vdef) {
    return MX::zz_substitute(ex, v, vdef);
  }

  /** \brief  Substitute variable v with expression vdef in an expression ex, preserving nodes */
  inline MX
  graph_substitute(const MX &ex, const std::vector<MX> &v, const std::vector<MX> &vdef) {
    return ex.zz_graph_substitute(v, vdef);
  }

  /** \brief  Substitute variable var with expression expr in
   * multiple expressions, preserving nodes  */
  inline std::vector<MX>
  graph_substitute(const std::vector<MX> &ex, const std::vector<MX> &v,
                   const std::vector<MX> &vdef) {
    return MX::zz_graph_substitute(ex, v, vdef);
  }

  /** \brief Inplace substitution
   * Substitute variables v out of the expressions vdef sequentially
   */
  inline void
  substituteInPlace(const std::vector<MX>& v,
                    std::vector<MX>& SWIG_INOUT(vdef), bool reverse=false) {
    return MX::zz_substituteInPlace(v, vdef, reverse);
  }

  /** \brief Inplace substitution with piggyback expressions
   * Substitute variables v out of the expressions vdef sequentially,
   * as well as out of a number of other expressions piggyback */
  inline void
  substituteInPlace(const std::vector<MX>& v,
                    std::vector<MX>& SWIG_INOUT(vdef),
                    std::vector<MX>& SWIG_INOUT(ex), bool reverse=false) {
    return MX::zz_substituteInPlace(v, vdef, ex, reverse);
  }

  /** \brief Extract shared subexpressions from an set of expressions */
  inline void extractShared(std::vector<MX>& ex,
                            std::vector<MX>& v, std::vector<MX>& vdef,
                            const std::string& v_prefix="v_",
                            const std::string& v_suffix="") {
    return MX::zz_extractShared(ex, v, vdef, v_prefix, v_suffix);
  }

  /** \brief Print compact, introducing new variables for shared subexpressions */
  inline void printCompact(const MX& ex, std::ostream &stream=CASADI_COUT) {
    return ex.zz_printCompact(stream);
  }

  ///@{
  /** \brief Calculate jacobian via source code transformation

      Uses casadi::MXFunction::jac
  */
  inline MX jacobian(const MX &ex, const MX &arg) { return ex.zz_jacobian(arg);}
  inline MX gradient(const MX &ex, const MX &arg) { return ex.zz_gradient(arg);}
  inline MX tangent(const MX &ex, const MX &arg) { return ex.zz_tangent(arg);}
  ///@}

  /** \brief Get all symbols contained in the supplied expression
  * Get all symbols on which the supplied expression depends
  * \see MXFunction::getFree()
  */
  inline std::vector<MX> getSymbols(const MX& e) { return e.zz_getSymbols();}

  /** \brief Get all symbols contained in the supplied expression
  * Get all symbols on which the supplied expression depends
  * \see MXFunction::getFree()
  */
  inline std::vector<MX> getSymbols(const std::vector<MX>& e) {
    return MX::zz_getSymbols(e);
  }

  /** \brief Check if expression depends on any of the arguments
    The arguments must be symbolic
  */
  inline bool dependsOn(const MX& ex, const std::vector<MX> &arg) { return ex.zz_dependsOn(arg);}

  /** \brief Expand MX graph to SXFunction call
  *
  *  Expand the given expression e, optionally
  *  supplying expressions contained in it at which expansion should stop.
  *
  */
  inline MX
  matrix_expand(const MX& e, const std::vector<MX> &boundary = std::vector<MX>()) {
    return MX::zz_matrix_expand(e, boundary);
  }

  /** \brief Expand MX graph to SXFunction call
  *
  *  Expand the given expression e, optionally
  *  supplying expressions contained in it at which expansion should stop.
  *
  */
  inline std::vector<MX>
  matrix_expand(const std::vector<MX>& e, const std::vector<MX> &boundary = std::vector<MX>()) {
    return MX::zz_matrix_expand(e, boundary);
  }

  /** \brief Solve a system of equations: A*x = b
  */
  inline MX
  solve(const MX& A, const MX& b, const std::string& lsolver = "symbolicqr",
        const Dictionary& dict = Dictionary()) { return A.zz_solve(b, lsolver, dict);}

  /** \brief Computes the Moore-Penrose pseudo-inverse
  *
  * If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity.
  * If the matrix A is slender (size2<size1), mul(pinv(A), A) is unity.
  *
  */
  inline MX pinv(const MX& A, const std::string& lsolver,
                 const Dictionary& dict = Dictionary()) { return A.zz_pinv(lsolver, dict);}

#endif // !defined(SWIG) || !defined(SWIGMATLAB)

/** @}
*/

} // namespace casadi

#endif // CASADI_MX_TOOLS_HPP
