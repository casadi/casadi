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

#ifndef SIMPLE_HOMOTOPY_NLP_SOLVER_HPP
#define SIMPLE_HOMOTOPY_NLP_SOLVER_HPP

#include "casadi/symbolic/function/homotopy_nlp_solver.hpp"

#include <casadi/nonlinear_programming/casadi_nonlinear_programming_export.h>

namespace casadi{

  class SimpleHomotopyNLPInternal;

  /**
     \brief Solving an NLP homotopy with regular NLP solvers

     \author Joris Gillis
     \date 2014
  */
  class CASADI_NONLINEAR_PROGRAMMING_EXPORT SimpleHomotopyNLPSolver : public HomotopyNLPSolver {
  public:
    /// Default constructor
    SimpleHomotopyNLPSolver();

    /// \brief Create an NLP solver instance
    explicit SimpleHomotopyNLPSolver(const Function& hnlp /**< nlp function: \f$ [\mathbb{R}^{n_x} \times \mathbb{R}^{n_p}] \mapsto [\mathbb{R} \times \mathbb{R}^{n_g}]\f$*/
                    );

    /// Access functions of the node
    SimpleHomotopyNLPInternal* operator->();
    const SimpleHomotopyNLPInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function
#ifdef SWIG
    %callback("%s_cb");
#endif
    static HomotopyNLPSolver creator(const Function& hnlp){ return SimpleHomotopyNLPSolver(hnlp);}
#ifdef SWIG
    %nocallback;
#endif


  };

} // namespace casadi

#endif //SIMPLE_HOMOTOPY_NLP_SOLVER_HPP
