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


#ifndef CASADI_DPLE_TO_LR_DPLE_HPP
#define CASADI_DPLE_TO_LR_DPLE_HPP

#include "../core/function/lr_dple_internal.hpp"
#include "../core/function/dple_internal.hpp"
#include <casadi/solvers/casadi_lrdplesolver_dple_export.h>

/** \defgroup plugin_LrDpleSolver_dple
 Solving the Low-Rank Discrete Lyapunov Equations with
 a Low-Rank Discrete Lyapunov Equations Solver

*/
/** \pluginsection{LrDpleSolver,dple} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{LrDpleSolver,dple}

   @copydoc LR_DLE_doc
   @copydoc plugin_LrDpleSolver_dple

       \author Joris Gillis
      \date 2014

  */
  class CASADI_LRDPLESOLVER_DPLE_EXPORT DpleToLrDple :
    public LrDpleInternal,
    public Adaptor<DpleToLrDple, DpleInternal>,
    public Wrapper<DpleToLrDple> {
  public:
    /** \brief  Constructor
     * \param st \structargument{LrDple}
     */
    DpleToLrDple(const LrDpleStructure & st,
                 const std::vector< std::vector<int> > &Hs);

    /** \brief  Destructor */
    virtual ~DpleToLrDple();

    /** \brief  Clone */
    virtual DpleToLrDple* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual DpleToLrDple* create(const LrDpleStructure& st,
                                 const std::vector< std::vector<int> > &Hs) const {
        return new DpleToLrDple(st, Hs);}

    /** \brief  Create a new DLE Solver */
    static LrDpleInternal* creator(const LrDpleStructure& st,
                                const std::vector< std::vector<int> > &Hs)
    { return new DpleToLrDple(st, Hs);}

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
    DpleSolver solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_DPLE_TO_LR_DPLE_HPP
