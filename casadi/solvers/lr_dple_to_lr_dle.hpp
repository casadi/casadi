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


#ifndef CASADI_LR_DPLE_TO_LR_DLE_HPP
#define CASADI_LR_DPLE_TO_LR_DLE_HPP

#include "../core/function/lr_dle_internal.hpp"
#include "../core/function/lr_dple_solver.hpp"
#include "../core/function/adaptor.hpp"
#include <casadi/solvers/casadi_lrdlesolver_lrdple_export.h>

/** \defgroup plugin_LrDleSolver_lrdple

 Solving the low-rank  Periodic Discrete Lyapunov Equations
 with a low-rank discrete Lyapunov solver

*/
/** \pluginsection{LrDleSolver,lrdple} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{LrDleSolver,lrdple}

   @copydoc DLE_doc
   @copydoc plugin_LrDleSolver_lrdple

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DLESOLVER_DPLE_EXPORT LrDpleToLrDle : public LrDleInternal,
    public Adaptor<LrDpleToLrDle, LrDleInternal> {
  public:
    /** \brief  Constructor
     *  \param[in] A  Sparsity of A
     *  \param[in] V  Sparsity of V
     */
    LrDpleToLrDle(const LrDleStructure& st, const std::vector<int> &Hs);

    /** \brief  Destructor */
    virtual ~LrDpleToLrDle();

    /** \brief  Clone */
    virtual LrDpleToLrDle* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual LrDpleToLrDle* create(const LrDleStructure& st, const std::vector<int> &Hs) const {
        return new LrDpleToLrDle(st, Hs);}

    /** \brief  Create a new DLE Solver */
    static DleInternal* creator(const LrDleStructure& st, const std::vector<int> &Hs)
    { return new LrDpleToLrDle(st, Hs);}

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

  };

} // namespace casadi
/// \endcond
#endif // CASADI_LR_DPLE_TO_LR_DLE_HPP
