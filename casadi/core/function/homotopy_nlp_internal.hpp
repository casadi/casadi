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


#ifndef CASADI_HOMOTOPY_NLP_SOLVER_INTERNAL_HPP
#define CASADI_HOMOTOPY_NLP_SOLVER_INTERNAL_HPP

#include "homotopy_nlp_solver.hpp"
#include "nlp_solver.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

/** \brief Homotopy NLP solver storage class
  \internal
  @copydoc HomotopyNlpSolver_doc
  \author Joris Gillis
  \date 2013-2014
*/
  class CASADI_CORE_EXPORT
  HomotopyNLPInternal : public FunctionInternal,
                        public PluginInterface<HomotopyNLPInternal> {
  public:
    /// Constructor
    HomotopyNLPInternal(const Function& hnlp);

    /// Destructor
    virtual ~HomotopyNLPInternal() = 0;

    /// Initialize
    virtual void init();

    /// The Homotopy NLP
    Function hnlp_;

    int nx_;
    int np_;
    int ng_;

    // Creator function for internal class
    typedef HomotopyNLPInternal* (*Creator)(const Function& hnlp);

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_NLP_SOLVER_INTERNAL_HPP
