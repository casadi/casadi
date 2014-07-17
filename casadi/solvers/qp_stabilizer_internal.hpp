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

#ifndef CASADI_QP_STABILIZER_INTERNAL_HPP
#define CASADI_QP_STABILIZER_INTERNAL_HPP

#include "casadi/core/function/stabilized_qp_solver_internal.hpp"
#include "casadi/core/function/stabilized_qp_solver.hpp"

#include <casadi/solvers/casadi_stabilizedqpsolver_qp_export.h>

/// \cond INTERNAL
namespace casadi {

  /** \brief Solved a stabilized QP using a standard QP solver

      \author Joris Gillis
      \date 2013
  */
  class CASADI_STABILIZEDQPSOLVER_QP_EXPORT QPStabilizerInternal
    : public StabilizedQpSolverInternal {
    friend class QPStabilizer;
  public:

    /** \brief Constructor */
    explicit QPStabilizerInternal(const std::vector<Sparsity> &st);

    /** \brief Destructor */
    virtual ~QPStabilizerInternal();

    /** \brief  Clone */
    virtual QPStabilizerInternal* clone() const { return new QPStabilizerInternal(*this);}

    /** \brief  Create a new Stabilized QP Solver */
    static StabilizedQpSolverInternal* creator(const QPStructure& st)
    { return new QPStabilizerInternal(st);}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief Initialize */
    virtual void init();

    /** \brief Solve the QP */
    virtual void evaluate();

    /** \brief Generate native code for debugging */
    virtual void generateNativeCode(std::ostream &file) const;

    /// Data members
    QpSolver qp_solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_QP_STABILIZER_INTERNAL_HPP

