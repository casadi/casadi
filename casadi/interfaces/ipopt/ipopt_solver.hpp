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

#ifndef IPOPT_SOLVER_HPP
#define IPOPT_SOLVER_HPP

#include "casadi/symbolic/function/nlp_solver.hpp"
#include <casadi/interfaces/ipopt/casadi_ipopt_interface_export.h>

namespace casadi{

  class IpoptInternal;

  // List from ipopt_internal.cpp
  /**
   * \brief interface to IPOPT NLP solver
   * @copydoc NLPSolver_doc
   *
   * When in warmstart mode, output NLP_SOLVER_LAM_X may be used as input
   *
   * NOTE: Even when max_iter == 0,  it is not guaranteed that
   *       input(NLP_SOLVER_X0) == output(NLP_SOLVER_X).
   *       Indeed if bounds on X or constraints are unmet, they will differ.
   *
   *  For a good tutorial on IPOPT, see 
   *  http://drops.dagstuhl.de/volltexte/2009/2089/pdf/09061.WaechterAndreas.Paper.2089.pdf
   *
   *  A good resource about the algorithms in IPOPT is: Wachter and L. T. Biegler,
   *  On the Implementation of an Interior-Point Filter Line-Search Algorithm for
   *  Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57,
   *  2006 (As Research Report RC 23149, IBM T. J. Watson Research Center, Yorktown, USA
   *
   * Caveats:
   *   * with default options, multipliers for the decision variables are wrong for equality
   *     constraints.
   *     Change the 'fixed_variable_treatment' to 'make_constraint' or 'relax_bounds' to obtain
   *     correct results.
   *
   *
   */
  class CASADI_IPOPT_INTERFACE_EXPORT IpoptSolver : public NLPSolver {
  public:
    /// Default constructor
    IpoptSolver();

    /// \brief Create an NLP solver instance
    explicit IpoptSolver(const Function& nlp
                         /**< nlp function: \f$ [\mathbb{R}^{n_x} \times \mathbb{R}^{n_p}] \mapsto
                          * [\mathbb{R} \times \mathbb{R}^{n_g}]\f$*/
                         );

    /** \brief Get the reduced Hessian.
     * Requires a patched sIPOPT installation, see CasADi documentation. */
    DMatrix getReducedHessian();

    /// Access functions of the node
    IpoptInternal* operator->();
    const IpoptInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function
#ifdef SWIG
    %callback("%s_cb");
#endif
    static NLPSolver creator(const Function& nlp){ return IpoptSolver(nlp);}
#ifdef SWIG
    %nocallback;
#endif

  };

} // namespace casadi

#endif //IPOPT_SOLVER_HPP
