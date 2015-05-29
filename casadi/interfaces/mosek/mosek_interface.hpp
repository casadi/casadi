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


#ifndef CASADI_MOSEK_INTERFACE_HPP
#define CASADI_MOSEK_INTERFACE_HPP

#include "casadi/core/function/socp_solver_internal.hpp"
#include <casadi/interfaces/mosek/casadi_socpsolver_mosek_export.h>

#include <mosek.h>

/** \defgroup plugin_SocpSolver_mosek
      Interface to the SOCP solver MOSEK
*/

/** \pluginsection{SocpSolver,mosek} */

/// \cond INTERNAL
namespace casadi {
  /** \brief \pluginbrief{SocpSolver,mosek}


   \author Niels van Duijkeren
   \date 2015
   *
   @copydoc SocpSolver_doc
   @copydoc plugin_SocpSolver_mosek
   * */
  class CASADI_SOCPSOLVER_MOSEK_EXPORT MosekInterface : public SocpSolverInternal {
  public:

    /** \brief Constructor */
    explicit MosekInterface(const std::vector<Sparsity> &st);

    /** \brief Clone */
    virtual MosekInterface* clone() const;

    /** \brief  Create a new SOCP Solver */
    static SocpSolverInternal* creator(const SOCPStructure& st)
    { return new MosekInterface(st);}

    /** \brief Destructor */
    virtual ~MosekInterface();

    /** \brief Initialize */
    virtual void init();

    /** \brief Solve the SOCP */
    virtual void evaluate();

    /** Get termination reason from flag */
    static const char* terminationReason(int flag);

    /** Get solution type from flag */
    static const char* solutionType(int flag);

    /// A documentation string
    static const std::string meta_doc;

  private:
    /** Generate dual of SOCP to obtain conic form required by Mosek */
    void convertToDualSocp();

    /** Get solution status from MOSEK solsta value */
    std::string solutionStatus(MSKsolstae& solsta);

    /** Get problem status from MOSEK prosta value */
    std::string problemStatus(MSKprostae& prosta);

    /** Dual problem data */
    std::vector<double>  dual_c_;
    std::vector<double>  dual_A_data_;
    std::vector<int>     dual_A_row_;
    std::vector<int>     dual_A_colind_;
    std::vector<double>  dual_b_;
    std::vector<double>  dual_yi_;
    std::vector<double>  dual_ti_;

    /** Indices of relevant bounds */
    std::vector<int>  primal_idx_lba_;     // Indices of lower bounded linear inequality constraints (LBA != -inf)
    std::vector<int>  primal_idx_uba_;     // Indices of upper bounded linear inequality constraints (UBA != inf)
    std::vector<int>  primal_idx_lbx_;     // Indices of simple lower bounds  (LBX != -inf)
    std::vector<int>  primal_idx_ubx_;     // Indices of simple upper bounds (UBX != inf)

    /** MOSEK variables */
    MSKenv_t   mosek_env_;
    MSKtask_t  mosek_task_;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_MOSEK_INTERFACE_HPP
