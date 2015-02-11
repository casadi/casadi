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


#ifndef CASADI_dple_to_lr_dple_HPP
#define CASADI_dple_to_lr_dple_HPP

#include "../core/function/dple_internal.hpp"
#include "../core/function/lr_dple_internal.hpp"
#include <casadi/solvers/casadi_dplesolver_lrdple_export.h>

/** \defgroup plugin_DpleSolver_lrdple
 Solving the Low-Rank Discrete Lyapunov Equations with
 a Low-Rank Discrete Lyapunov Equations Solver

*/
/** \pluginsection{DpleSolver,lrdple} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{DpleSolver,lrdple}

   @copydoc DPLE_doc
   @copydoc plugin_DpleSolver_lrdple

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLESOLVER_LRDPLE_EXPORT DpleToLrDple : public DpleInternal,
    public Adaptor<DpleToLrDple, LrDpleInternal>,
    public Wrapper<DpleToLrDple> {
  public:
    /** \brief  Constructor
     * \param st \structargument{LrDple}
     */
    DpleToLrDple(const DpleStructure & st);

    /** \brief  Destructor */
    virtual ~DpleToLrDple();

    /** \brief  Clone */
    virtual DpleToLrDple* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual DpleToLrDple* create(const DpleStructure& st) const {
        return new DpleToLrDple(st);}

    /** \brief  Create a new DLE Solver */
    static DpleInternal* creator(const DpleStructure& st)
    { return new DpleToLrDple(st);}

    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const {}

    /** \brief  evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    /** Generate a function that calculates \a nfwd forward derivatives */
    virtual Function getDerivativeFwd(int nfwd);

    /** Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function getDerivativeAdj(int nadj);

    /// A documentation string
    static const std::string meta_doc;

    /// Solve with
    LrDpleSolver solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_dple_to_lr_dple_HPP
