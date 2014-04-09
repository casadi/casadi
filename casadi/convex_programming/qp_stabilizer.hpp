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

#ifndef QP_STABILIZER_HPP
#define QP_STABILIZER_HPP

#include "casadi/symbolic/function/stabilized_qp_solver.hpp"

#include <casadi/convex_programming/casadi_convex_programming_export.h>

namespace casadi {


// Forward declaration of internal class
class QPStabilizerInternal;

  /** \brief IPOPT QP Solver for quadratic programming

   @copydoc StabilizedQPSolver_doc

   \author Joris Gillis
   \date 2013
  */
class CASADI_CONVEX_PROGRAMMING_EXPORT QPStabilizer : public StabilizedQPSolver {
public:

  /** \brief  Default constructor */
  QPStabilizer();


  /** \brief Constructor
  *  \param st Problem structure
  *  \copydoc scheme_QPStruct
  */
  explicit QPStabilizer(const QPStructure & st);

  /** \brief  Access functions of the node */
  QPStabilizerInternal* operator->();
  const QPStabilizerInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static StabilizedQPSolver creator(const QPStructure & st){ return QPStabilizer(st);}
  #ifdef SWIG
  %nocallback;
  #endif

  /// Access underlying QP solver
  QPSolver & getSolver();

};


} // namespace casadi

#endif //QP_STABILIZER_HPP

