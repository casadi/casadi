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


#ifndef CASADI_CONDENSING_INDEF_DPLE_INTERNAL_HPP
#define CASADI_CONDENSING_INDEF_DPLE_INTERNAL_HPP

#include "../core/function/dple_internal.hpp"
#include "../core/function/dle_internal.hpp"
#include <casadi/solvers/casadi_dplesolver_condensing_export.h>

/** \defgroup plugin_DpleSolver_condensing
 Solving the Discrete Periodic Lyapunov Equations by
 condensing the entire period to a single Discrete Lyapunov Equation

*/
/** \pluginsection{DpleSolver,condensing} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{DpleSolver,condensing}

    The condensing approach transforms the discrete periodic Lyapunov equations
    into one small and dense discrete Lyapunov equation.

    An example for a period K=3  is given below.

    The discrete periodic Lyapunov equations:

    \verbatim
      P1 = A0 P0 A0^T + V0
      P2 = A1 P1 A1^T + V1
      P0 = A2 P2 A2^T + V2
    \endverbatim

    Substitute the equations to obtain:

    \verbatim
    P0 = A2 ( A1 ( A0 P0 A0^T + V0 ) A1^T + V1 ) A2^T + V2
    \endverbatim

    Hence, we have
    \verbatim
    P = A P A^T +R

     A = A2 A1 A0
     R = A2 A1 V0 A1^T A2^T + A2 V1 A2^T + V2
    \endverbatim

    The desired outputs are reconstructed as:
    \verbatim
    P0 = P
    P1 = A0 P0 A0^T + V0
    P2 = A1 P1 A1^T + V1
    \endverbatim

   @copydoc DPLE_doc
   @copydoc plugin_DpleSolver_condensing

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLESOLVER_CONDENSING_EXPORT CondensingIndefDpleInternal : public DpleInternal,
    public Adaptor<CondensingIndefDpleInternal, DleInternal>,
    public Wrapper<CondensingIndefDpleInternal> {
  public:
    /** \brief  Constructor
     * \param st \structargument{Dple}
     */
    CondensingIndefDpleInternal(const DpleStructure& st);

    /** \brief  Destructor */
    virtual ~CondensingIndefDpleInternal();

    /** \brief  Clone */
    virtual CondensingIndefDpleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual CondensingIndefDpleInternal* create(const DpleStructure& st) const {
        return new CondensingIndefDpleInternal(st);}

    /** \brief  Create a new DPLE Solver */
    static DpleInternal* creator(const DpleStructure& st) {
        return new CondensingIndefDpleInternal(st);}

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

  private:

    /// State space dimension
    int n_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_CONDENSING_INDEF_DPLE_INTERNAL_HPP
