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


#ifndef CASADI_FIXED_SMITH_DLE_INTERNAL_HPP
#define CASADI_FIXED_SMITH_DLE_INTERNAL_HPP

#include "../core/function/dle_internal.hpp"
#include <casadi/solvers/casadi_dlesolver_fixed_smith_export.h>

/** \defgroup plugin_DleSolver_fixed_smith
 Solving the Discrete Lyapunov Equations
 with a fixed number of smith iterations.
 
 @copdyoc DleSolversmith

*/
/** \defgroup DleSolversmith

 This plugin uses Smith iterations.
 
 Th basic idea is to exploit the fact that the discrete 
 algebraic Lyapunov operator
 f(X) = AXA^T + V
 has a fixed point when A is stable.
 
 The pure Smith iterations are:
 
 \verbatim
 
 X_{-1} = 0
 X_0 = V
 k = 0
 while ||X_k − X_{k−1} || < do
   X_{k+1} = A X_k A^T + V
   k += 1
 end
 
 P = X_k
 \endverbatim
 
 With frequency doubling, we have:
 
 \verbatim
 
 X_{-1} = 0
 X_0 = V
 V_0 = V
 A_0 = A
 k = 0
 while ||X_k − X_{k−1} || < do
   X_{k+1} = A_k X_k A_k^T + V_k
   V_{k+1} = A_k V_k A_k^T + V_k
   A_{k+1} = A_k A_k
   k += 1
 end
 
 P = X_k
 \endverbatim


*/
/** \pluginsection{DleSolver,fixed_smith} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{DleSolver,fixed_smith}

   @copydoc DLE_doc
   @copydoc plugin_DleSolver_fixed_smith

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DLESOLVER_FIXED_SMITH_EXPORT FixedSmithDleInternal : public DleInternal,
    public Wrapper<FixedSmithDleInternal> {
  public:
    /** \brief  Constructor
     * \param st \structargument{Dle}
     * \param Hs Column-sizes of H_i
     */
    FixedSmithDleInternal(const DleStructure& st);

    /** \brief  Destructor */
    virtual ~FixedSmithDleInternal();

    /** \brief  Clone */
    virtual FixedSmithDleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual FixedSmithDleInternal* create(const DleStructure& st) const {
        return new FixedSmithDleInternal(st);}

    /** \brief  Create a new DLE Solver */
    static DleInternal* creator(const DleStructure& st)
    { return new FixedSmithDleInternal(st);}

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

    /// Number of Smith iterations
    int iter_;

    /// Frequency doubling?
    bool freq_doubling_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_FIXED_SMITH_DLE_INTERNAL_HPP
