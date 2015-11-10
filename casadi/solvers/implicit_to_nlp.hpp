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

#include "casadi/core/function/nlsol.hpp"
#include "casadi/core/function/nlpsol.hpp"

#include <casadi/solvers/casadi_nlsol_nlpsol_export.h>

/** \defgroup plugin_Nlsol_nlp
  Use an Nlpsol as Nlsol solver
*/
/** \pluginsection{Nlsol,nlp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief  \pluginbrief{Nlsol,nlp}

   @copydoc Nlsol_doc
   @copydoc plugin_Nlsol_nlp

   \author Joris Gillis
   \date 2012
  */
  class CASADI_NLSOL_NLPSOL_EXPORT ImplicitToNlp : public Nlsol {
  public:
    /** \brief  Constructor */
    explicit ImplicitToNlp(const std::string& name, const Function& f);

    /** \brief  Destructor */
    virtual ~ImplicitToNlp();

    /** \brief  Create a new Nlsol */
    static Nlsol* creator(const std::string& name, const Function& f) {
      return new ImplicitToNlp(name, f);
    }

    // Get name of the plugin
    virtual const char* plugin_name() const { return "nlpsol";}

    /** \brief  Initialize */
    virtual void init();

    /// Solve the system of equations and calculate derivatives
    virtual void eval(const double** arg, double** res, int* iw, double* w, void* mem);

    /// A documentation string
    static const std::string meta_doc;

    /// NLP solver
    Function solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_IMPLICIT_TO_NLP_HPP
