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

#ifndef OOQP_INTERNAL_HPP
#define OOQP_INTERNAL_HPP

#include "casadi/symbolic/function/qp_solver_internal.hpp"
#include "ooqp_solver.hpp"

/// \cond INTERNAL
namespace casadi{

  /** \brief Internal class for OOQPSolver
   *
   @copydoc QPSolver_doc
   * */
  class CASADI_OOQP_INTERFACE_EXPORT OOQPInternal : public QPSolverInternal {
    friend class OOQPSolver;
  public:

    /** \brief  Constructor */
    explicit OOQPInternal();

    /** \brief  Clone */
    virtual OOQPInternal* clone() const{ return new OOQPInternal(*this);}

    /** \brief  Create a new Solver */
    explicit OOQPInternal(const std::vector<Sparsity>& st);

    /** \brief  Destructor */
    virtual ~OOQPInternal();

    /** \brief  Initialize */
    virtual void init();

    /// Solve the QP
    virtual void evaluate();

    /// Throw error
    static const char* errFlag(int flag);

    /// Print an OOQP bounds vector
    static std::string printBounds(const std::vector<double>& b, const std::vector<char>& ib, int n, const char *sign);

    /// Problem data (vectors)
    std::vector<double> c_, bA_, xlow_, xupp_, clow_, cupp_, x_, gamma_, phi_, y_, z_, lambda_, pi_;

    /// Type of bounds
    std::vector<char> ixlow_, ixupp_, iclow_, icupp_;

    /// Problem data (matrices)
    std::vector<double> dQ_, dA_, dC_;

    /// Sparsities of matrices
    std::vector<int> irowQ_, jcolQ_, irowA_, jcolA_, irowC_, jcolC_;

    // Variable/constraint index
    std::vector<int> x_index_, c_index_;

    // Parameters
    std::vector<double> p_;

    // Transpose of linear constraints
    DMatrix AT_;
    std::vector<int> AT_tmp_;

    // Print level
    int print_level_;

    // Tolerances
    double mutol_, artol_;
  };

} // namespace casadi

/// \endcond
#endif //OOQP_INTERNAL_HPP

