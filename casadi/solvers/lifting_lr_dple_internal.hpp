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


#ifndef CASADI_LIFTING_LR_DPLE_INTERNAL_HPP
#define CASADI_LIFTING_LR_DPLE_INTERNAL_HPP

#include "../core/function/lr_dple_internal.hpp"
#include "../core/function/lr_dle_internal.hpp"
#include <casadi/solvers/casadi_lrdplesolver_lifting_export.h>

/** \defgroup plugin_LrDpleSolver_lifting
 Solving the Discrete Periodic Lyapunov Equations by
 lifting the entire period to a single Discrete Lyapunov Equation

*/
/** \pluginsection{LrDpleSolver,lifting} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{LrDpleSolver,lifting}

   @copydoc DPLE_doc
   @copydoc plugin_LrDpleSolver_lifting


    The lifting approach transforms the discrete periodic Lyapunov equations
    into one large and sparse discrete Lyapunov equation.

    An example for a period K=3  is given below.

    The discrete periodic Lyapunov equations:

    \verbatim
      P1 = A0 P0 A0^T + V0
      P2 = A1 P1 A1^T + V1
      P0 = A2 P2 A2^T + V2
    \endverbatim

    The lifted form (type A):
    \verbatim
                                                      T
     | P0     |   |     A2 |   | P0     |   |     A2 |   | V2     |
     |   P1   | = | A0     | . |   P1   | . | A0     | + |   V0   |
     |     P2 |   |   A1   |   |     P2 | . |   A1   |   |     V1 |
    \endverbatim

    The lifted form (type B):
    \verbatim
                                                      T
     | P2     |   |   A1   |   | P2     |   |   A1   |   | V1     |
     |   P1   | = |     A0 | . |   P1   | . |     A0 | + |   V0   |
     |     P0 |   | A2     |   |     P0 | . | A2     |   |     V2 |
    \endverbatim

       \author Joris Gillis
      \date 2014

  */
  class CASADI_LRDPLESOLVER_LIFTING_EXPORT LiftingLrDpleInternal : public LrDpleInternal,
    public Adaptor<LiftingLrDpleInternal, LrDleInternal>,
    public Wrapper<LiftingLrDpleInternal> {
  public:
    /** \brief  Constructor
     * \param st \structargument{LrDple}
     */
    LiftingLrDpleInternal(const LrDpleStructure & st,
                             const std::vector< std::vector<int> > &Hs);

    /** \brief  Destructor */
    virtual ~LiftingLrDpleInternal();

    /** \brief  Clone */
    virtual LiftingLrDpleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual LiftingLrDpleInternal* create(const LrDpleStructure & st,
      const std::vector< std::vector<int> > &Hs) const {
        return new LiftingLrDpleInternal(st, Hs);}

    /** \brief  Create a new DPLE Solver */
    static LrDpleInternal* creator(const LrDpleStructure & st,
      const std::vector< std::vector<int> > &Hs) {
        return new LiftingLrDpleInternal(st, Hs);}

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
    LrDleSolver solver_;

  private:

    /// State space dimension
    int n_;

    /// Form of the lifted equations
    int form_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_LIFTING_INDEF_DPLE_INTERNAL_HPP
