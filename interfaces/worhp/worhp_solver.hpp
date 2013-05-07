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

#ifndef WORHP_SOLVER_HPP
#define WORHP_SOLVER_HPP

#include "symbolic/fx/nlp_solver.hpp"

namespace CasADi{
  
  class WorhpInternal;
  
  // List from worhp_internal.cpp
  /**
   *
   * \brief interface to WORHP NLP solver
   * @copydoc NLPSolver_doc

   */
  class WorhpSolver : public NLPSolver {
  public:
    /// Default constructor
    WorhpSolver();

    /// \brief Create an NLP solver instance (legacy syntax)
    explicit WorhpSolver(const FX& F, /**< objective function: \f$ [\mathbb{R}^{n_x}] \mapsto [\mathbb{R}]\f$*/
                         const FX& G  /**< constraint function \f$ [\mathbb{R}^{n_x}] \mapsto [\mathbb{R}^{n_g}]\f$ */
                         );

    /// \brief Create an NLP solver instance
    explicit WorhpSolver(const FX& nlp /**< nlp function: \f$ [\mathbb{R}^{n_x} \times \mathbb{R}^{n_p}] \mapsto [\mathbb{R} \times \mathbb{R}^{n_g}]\f$*/
                         );

    /// Access functions of the node
    WorhpInternal* operator->();
    const WorhpInternal* operator->() const;

    /// Read options from worhp parameter xml
    void setOptionsFromFile(const std::string & file);

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function 
#ifdef SWIG
    %callback("%s_cb");
#endif
    static NLPSolver creator(const FX& nlp){ return WorhpSolver(nlp);}
#ifdef SWIG
    %nocallback;
#endif
  };

} // namespace CasADi

#endif //WORHP_SOLVER_HPP
