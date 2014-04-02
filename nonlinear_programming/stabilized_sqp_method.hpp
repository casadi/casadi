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

#ifndef STABILIZED_SQP_METHOD_HPP
#define STABILIZED_SQP_METHOD_HPP

#include "symbolic/function/nlp_solver.hpp"

namespace CasADi{
  
  class StabilizedSQPInternal;
  
  /**
     \brief Stabilized Sequential Quadratic Programming method.

     \author Slava Kung
     \date 2013
  */
  class StabilizedSQPMethod : public NLPSolver {
  public:
    /// Default constructor
    StabilizedSQPMethod();

#ifndef WITHOUT_PRE_1_9_X
    /// \brief [DEPRECATED] Create an NLP solver instance (legacy syntax)
    explicit StabilizedSQPMethod(const Function& F, /**< objective function: \f$ [\mathbb{R}^{n_x}] \mapsto [\mathbb{R}]\f$*/
                       const Function& G  /**< constraint function \f$ [\mathbb{R}^{n_x}] \mapsto [\mathbb{R}^{n_g}]\f$ */
                    );
#endif

    /// \brief Create an NLP solver instance
    explicit StabilizedSQPMethod(const Function& nlp /**< nlp function: \f$ [\mathbb{R}^{n_x} \times \mathbb{R}^{n_p}] \mapsto [\mathbb{R} \times \mathbb{R}^{n_g}]\f$*/
                    );

    /// Access functions of the node
    StabilizedSQPInternal* operator->();
    const StabilizedSQPInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function 
#ifdef SWIG
    %callback("%s_cb");
#endif
    static NLPSolver creator(const Function& nlp){ return StabilizedSQPMethod(nlp);}
#ifdef SWIG
    %nocallback;
#endif

    /// Access the StabilizedQPSolver used internally
    const StabilizedQPSolver getStabilizedQPSolver() const;
    
  };

} // namespace CasADi

#endif //STABILIZED_SQP_METHOD_HPP
