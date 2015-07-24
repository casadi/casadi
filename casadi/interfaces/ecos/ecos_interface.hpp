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


#ifndef CASADI_ECOS_INTERFACE_HPP
#define CASADI_ECOS_INTERFACE_HPP

#include "casadi/core/function/socp_solver_internal.hpp"
#include <casadi/interfaces/ecos/casadi_socpsolver_ecos_export.h>

#include <ecos.h>

/** \defgroup plugin_SocpSolver_ecos
      Interface to the SOCP solver ECOS
*/

/** \pluginsection{SocpSolver,ecos} */

/// \cond INTERNAL
namespace casadi {
  /** \brief \pluginbrief{SocpSolver,ecos}


   \author Niels van Duijkeren
   \date 2015
   *
   @copydoc SocpSolver_doc
   @copydoc plugin_SocpSolver_ecos
   * */
  class CASADI_SOCPSOLVER_ECOS_EXPORT EcosInterface : public SocpSolverInternal {
  public:

    /** \brief Constructor */
    explicit EcosInterface(const std::vector<Sparsity> &st);

    /** \brief Clone */
    virtual EcosInterface* clone() const;

    /** \brief  Create a new SOCP Solver */
    static SocpSolverInternal* creator(const SOCPStructure& st)
    { return new EcosInterface(st);}

    /** \brief Destructor */
    virtual ~EcosInterface();

    /** \brief Initialize */
    virtual void init();

    /** \brief Solve the SOCP */
    virtual void evaluate();

    /// A documentation string
    static const std::string meta_doc;

  private:

    /// Static ECOS input: G
    std::vector<pfloat> ecos_Gpr_vec_;
    std::vector<idxint> ecos_Gjc_vec_;
    std::vector<idxint> ecos_Gir_vec_;

    /// Static ECOS input: h
    std::vector<pfloat> ecos_h_vec_;

    /// Static ECOS input: q
    std::vector<idxint> ecos_q_vec_;

    /// ECOS input: c
    std::vector<pfloat> ecos_c_vec_;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_ECOS_INTERFACE_HPP
