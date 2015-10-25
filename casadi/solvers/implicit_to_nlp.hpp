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


#ifndef CASADI_IMPLICIT_TO_NLP_HPP
#define CASADI_IMPLICIT_TO_NLP_HPP

#include "casadi/core/function/rootfinder.hpp"
#include "casadi/core/function/nlp_solver.hpp"
#include "casadi/core/function/linear_solver.hpp"

#include <casadi/solvers/casadi_rootfinder_nlp_export.h>

/** \defgroup plugin_Rootfinder_nlp
  Use an NlpSolver as Rootfinder solver
*/
/** \pluginsection{Rootfinder,nlp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief  \pluginbrief{Rootfinder,nlp}

   @copydoc Rootfinder_doc
   @copydoc plugin_Rootfinder_nlp

   \author Joris Gillis
   \date 2012
  */
  class CASADI_ROOTFINDER_NLP_EXPORT QpToImplicit : public Rootfinder,
    public Adaptor<QpToImplicit, NlpSolver> {
  public:
    /** \brief  Constructor */
    explicit QpToImplicit(const std::string& name, const Function& f);

    /** \brief  Destructor */
    virtual ~QpToImplicit();

    /** \brief  Create a new Rootfinder */
    virtual QpToImplicit* create(const std::string& name, const Function& f) const {
      return new QpToImplicit(name, f);
    }

    /** \brief  Create a new Rootfinder */
    static Rootfinder* creator(const std::string& name, const Function& f) {
      return new QpToImplicit(name, f);
    }

    // Get name of the plugin
    virtual const char* plugin_name() const { return "nlp";}

    /** \brief  Initialize */
    virtual void init();

    /** \brief  Solve the nonlinear system of equations */
    virtual void solveNonLinear();

    /// A documentation string
    static const std::string meta_doc;

    /// Solve with
    Function solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_IMPLICIT_TO_NLP_HPP
