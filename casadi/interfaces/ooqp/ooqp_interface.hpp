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


#ifndef CASADI_OOQP_INTERFACE_HPP
#define CASADI_OOQP_INTERFACE_HPP

#include "casadi/core/function/qp_solver_internal.hpp"
#include <casadi/interfaces/ooqp/casadi_qpsolver_ooqp_export.h>

/** \defgroup plugin_QpSolver_ooqp
 Interface to the OOQP Solver for quadratic programming
  The current implementation assumes that OOQP is configured with the MA27 sparse linear solver.

  NOTE: when doing multiple calls to evaluate(), check if you need to reInit();
*/

/** \pluginsection{QpSolver,ooqp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{QpSolver,ooqp}

      @copydoc QpSolver_doc
      @copydoc plugin_QpSolver_ooqp

  */
  class CASADI_QPSOLVER_OOQP_EXPORT OoqpInterface : public QpSolverInternal {
  public:

    /** \brief  Constructor */
    explicit OoqpInterface();

    /** \brief  Clone */
    virtual OoqpInterface* clone() const { return new OoqpInterface(*this);}

    /** \brief  Create a new Solver */
    explicit OoqpInterface(const std::vector<Sparsity>& st);

    /** \brief  Create a new QP Solver */
    static QpSolverInternal* creator(const QPStructure& st)
    { return new OoqpInterface(st);}

    /** \brief  Destructor */
    virtual ~OoqpInterface();

    /** \brief  Initialize */
    virtual void init();

    /// Solve the QP
    virtual void evaluate();

    /// Throw error
    static const char* errFlag(int flag);

    /// Print an OOQP bounds vector
    static std::string printBounds(const std::vector<double>& b,
                                   const std::vector<char>& ib, int n, const char *sign);

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

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_OOQP_INTERFACE_HPP

