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


#ifndef CASADI_QPOASES_INTERFACE_HPP
#define CASADI_QPOASES_INTERFACE_HPP

#include "casadi/core/function/qp_solver_internal.hpp"
#include <casadi/interfaces/qpoases/casadi_qpsolver_qpoases_export.h>
#include <qpOASES.hpp>

/** \defgroup plugin_QpSolver_qpoases
Interface to QPOases Solver for quadratic programming

*/

/** \pluginsection{QpSolver,qpoases} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{QpSolver,qpoases}
   *
   * @copydoc QPSolver_doc
   * @copydoc plugin_QpSolver_qpoases
   *
   * \author Joris Gillis, Joel Andersson
   * \date 2011
   *
   * */
class CASADI_QPSOLVER_QPOASES_EXPORT QpoasesInterface : public QpSolverInternal {
public:
  /** \brief  Constructor */
  explicit QpoasesInterface();

  /** \brief  Clone */
  virtual QpoasesInterface* clone() const;

  /** \brief  Create a new QP Solver */
  static QpSolverInternal* creator(const QPStructure& st)
  { return new QpoasesInterface(st);}

  /** \brief  Create a new Solver */
  explicit QpoasesInterface(const std::vector<Sparsity>& st);

  /** \brief  Destructor */
  virtual ~QpoasesInterface();

  /** \brief  Initialize */
  virtual void init();

  virtual void evaluate();

  /// A documentation string
  static const std::string meta_doc;

  protected:

    /// QP Solver
    union {
      qpOASES::QProblemB *qp_;
    };

    ///@{
      /// Convert between qpOASES types and standard types
    static bool BooleanType_to_bool(qpOASES::BooleanType b);
    static qpOASES::BooleanType bool_to_BooleanType(bool b);
    static std::string SubjectToStatus_to_string(qpOASES::SubjectToStatus b);
    static qpOASES::SubjectToStatus string_to_SubjectToStatus(std::string b);
    static std::string PrintLevel_to_string(qpOASES::PrintLevel b);
    static qpOASES::PrintLevel string_to_PrintLevel(std::string b);
    ///@}

    /// Number of working set recalculations
    int max_nWSR_;

    /// CPUtime for initialization
    double max_cputime_;

    /// Throw error
    static void qpoases_error(const std::string& module, int flag);

    /// Has qpOASES been called once?
    bool called_once_;

    /// Dense data for H and A
    std::vector<double> h_data_;
    std::vector<double> a_data_;

    /// Temporary vector holding the dual solution
    std::vector<double> dual_;

    /// Get qpOASES error message
    static std::string getErrorMessage(int flag);

};

} // namespace casadi

/// \endcond
#endif // CASADI_QPOASES_INTERFACE_HPP
