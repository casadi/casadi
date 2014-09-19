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

#ifndef CASADI_DPLE_TO_DLE_HPP
#define CASADI_DPLE_TO_DLE_HPP

#include "../core/function/dle_internal.hpp"
#include "../core/function/dple_solver.hpp"
#include <casadi/solvers/casadi_dlesolver_dple_export.h>

/** \defgroup plugin_DleSolver_dple
 Solving the Discrete Lyapunov Equations with Periodic Solver

*/
/** \pluginsection{DleSolver,dple} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{DleSolver,dple}

   @copydoc DLE_doc
   @copydoc plugin_DleSolver_dple

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DLESOLVER_DPLE_EXPORT DpleToDle : public DleInternal {
  public:
    /** \brief  Constructor
     *  \param[in] A  Sparsity of A
     *  \param[in] V  Sparsity of V
     */
    DpleToDle(const Sparsity & A, const Sparsity &V);

    /** \brief  Destructor */
    virtual ~DpleToDle();

    /** \brief  Clone */
    virtual DpleToDle* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual DpleToDle* create(const Sparsity & A,
                                            const Sparsity &V) const {
        return new DpleToDle(A, V);}

    /** \brief  Create a new DLE Solver */
    static DleInternal* creator(const Sparsity & A,
                                 const Sparsity & V)
    { return new DpleToDle(A, V);}

    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const {}

    /** \brief  evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    /** \brief Generate a function that calculates \a nfwd forward derivatives
     and \a nadj adjoint derivatives
    */
    virtual Function getDerivative(int nfwd, int nadj);

    /// A documentation string
    static const std::string meta_doc;

  private:
    
    /// The solver used internally 
    DpleSolver dplesolver_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_DPLE_TO_DLE_HPP
