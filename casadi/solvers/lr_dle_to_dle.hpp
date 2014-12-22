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


#ifndef CASADI_LR_DLE_TO_DLE_HPP
#define CASADI_LR_DLE_TO_DLE_HPP

#include "../core/function/lr_dle_internal.hpp"
#include "../core/function/dle_internal.hpp"
#include <casadi/solvers/casadi_lrdle_dle_export.h>

/** \defgroup plugin_LrDleSolver_dle
 Solving the Low-Rank Discrete Lyapunov Equations with a
 regular Discrete Lyapunov Equations solver

*/
/** \pluginsection{LrDleSolver,dle} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{LrDleSolver,dle}

   @copydoc LR_DLE_doc
   @copydoc plugin_LrDleSolver_dle

       \author Joris Gillis
      \date 2014

  */
  class CASADI_LRDLE_DLE_EXPORT LrDleToDle : public LrDleInternal,
    public Adaptor<LrDleToDle, DleInternal>,
    public Wrapper<LrDleToDle> {
  public:
    /** \brief  Constructor
     * \param st \structargument{LrDle}
     * \param Hs Column-sizes of H_i
     */
    LrDleToDle(const LrDleStructure& st);

    /** \brief  Destructor */
    virtual ~LrDleToDle();

    /** \brief  Clone */
    virtual LrDleToDle* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual LrDleToDle* create(const LrDleStructure& st) const {
        return new LrDleToDle(st);}

    /** \brief  Create a new DLE Solver */
    static LrDleInternal* creator(const LrDleStructure& st)
    { return new LrDleToDle(st);}

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

    /// Solve with
    DleSolver solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_LR_DLE_TO_DLE_HPP
