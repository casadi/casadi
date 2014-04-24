/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef SOCP_QCQP_SOLVER_HPP
#define SOCP_QCQP_SOLVER_HPP

#include "casadi/core/function/qcqp_solver.hpp"

#include <casadi/convex_programming/casadi_convex_programming_export.h>

namespace casadi {


// Forward declaration of internal class
class SOCPQCQPInternal;

  /** \brief SOCP QCQP Solver for quadratic programming
   *
   *  Note: this implementation relies on Cholesky decomposition:
   *        <tt>Chol(H) = L ->  H = LL'</tt> with L lower triangular
   *   This requires Pi, H to be positive definite. Positive semi-definite is not sufficient.
   *    Notably, H==0  will not work.
   *
   *  A better implementation would rely on matrix square root, but we need
   *  singular value decomposition to implement that.
   *
   *
   * This implementation makes use of the epigraph reformulation:
   * \verbatim 
   *  min f(x)
   *    x
   *
   *   min  t
   *    x, t  f(x) <= t
   * \endverbatim
   *
   *  This implementation makes use of the following identity:
   * \verbatim
   *  || Gx+h||_2 <= e'x + f
   *
   *  x'(G'G - ee')x + (2 h'G - 2 f e') x + h'h - f <= 0
   * \endverbatim
   *    where we put e = [0 0 ... 1]  for the quadratic constraint
   *    arising from the epigraph reformulation and e==0 for all other 
   *    quadratic constraints.

   @copydoc QCQPSolver_doc

   \author Joris Gillis
   \date 2013
  */
class CASADI_CONVEX_PROGRAMMING_EXPORT SOCPQCQPSolver : public QCQPSolver {
public:

  /** \brief  Default constructor */
  SOCPQCQPSolver();


  /** \brief Constructor
  *  \param st Problem structure
  *  \copydoc scheme_QCQPStruct
  */
  explicit SOCPQCQPSolver(const QCQPStructure & st);

  /** \brief  Access functions of the node */
  SOCPQCQPInternal* operator->();
  const SOCPQCQPInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static QCQPSolver creator(const QCQPStructure & st) { return SOCPQCQPSolver(st);}
  #ifdef SWIG
  %nocallback;
  #endif

  /// Access underlying SOCP solver
  SOCPSolver & getSolver();

};


} // namespace casadi

#endif //SOCP_QCQP_SOLVER_HPP
