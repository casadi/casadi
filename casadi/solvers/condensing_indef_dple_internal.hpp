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

#ifndef CASADI_CONDENSING_INDEF_DPLE_INTERNAL_HPP
#define CASADI_CONDENSING_INDEF_DPLE_INTERNAL_HPP

#include "../core/function/dple_internal.hpp"
#include "../core/function/dle_solver.hpp"
#include <casadi/solvers/casadi_dplesolver_condensing_export.h>

/** \defgroup plugin_DpleSolver_condensing
 Solving the Discrete Periodic Lyapunov Equations by
 condensing the entire period to a single Discrete Lyapunov Equation

*/
/** \pluginsection{DpleSolver,condensing} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{DpleSolver,condensing}

   @copydoc DPLE_doc
   @copydoc plugin_DpleSolver_condensing

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLESOLVER_CONDENSING_EXPORT CondensingIndefDpleInternal : public DpleInternal {
  public:
    /** \brief  Constructor
     *  \param[in] A  List of sparsities of A_i
     *  \param[in] V  List of sparsities of V_i
     */
    CondensingIndefDpleInternal(const std::vector< Sparsity > & A, const std::vector< Sparsity > &V);

    /** \brief  Destructor */
    virtual ~CondensingIndefDpleInternal();

    /** \brief  Clone */
    virtual CondensingIndefDpleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual CondensingIndefDpleInternal* create(const std::vector< Sparsity > & A,
                                            const std::vector< Sparsity > &V) const {
        return new CondensingIndefDpleInternal(A, V);}

    /** \brief  Create a new DPLE Solver */
    static DpleInternal* creator(const std::vector< Sparsity >& A,
                                 const std::vector< Sparsity >& V)
    { return new CondensingIndefDpleInternal(A, V);}

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
    /// Main implementation as MXFunction
    Function f_;
    
    /// DleSolver solving the condensed Lyapunov form
    DleSolver dlesolver_;
    
    /// State space dimension
    int n_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_CONDENSING_INDEF_DPLE_INTERNAL_HPP
