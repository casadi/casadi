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


#ifndef CASADI_DLE_TO_LR_DLE_HPP
#define CASADI_DLE_TO_LR_DLE_HPP

#include "../core/function/lr_dle_internal.hpp"
#include "../core/function/dle_internal.hpp"
#include <casadi/solvers/casadi_lrdlesolver_dle_export.h>

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
  class CASADI_LRDLESOLVER_DLE_EXPORT DleToLrDle : public LrDleInternal,
    public Adaptor<DleToLrDle, DleInternal>,
    public Wrapper<DleToLrDle> {
  public:
    /** \brief  Constructor
     * \param st \structargument{LrDle}
     * \param Hs Column-sizes of H_i
     */
    DleToLrDle(const LrDleStructure& st, const std::vector<int> &Hs);

    /** \brief  Destructor */
    virtual ~DleToLrDle();

    /** \brief  Clone */
    virtual DleToLrDle* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual DleToLrDle* create(const LrDleStructure& st, const std::vector<int> &Hs) const {
        return new DleToLrDle(st, Hs);}

    /** \brief  Create a new DLE Solver */
    static LrDleInternal* creator(const LrDleStructure& st, const std::vector<int> &Hs)
    { return new DleToLrDle(st, Hs);}

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
#endif // CASADI_DLE_TO_LR_DLE_HPP
