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

#ifndef DPLE_INTERNAL_HPP
#define DPLE_INTERNAL_HPP

#include "dple_solver.hpp"
#include "dple_internal.hpp"
#include "../symbolic/function/function_internal.hpp"

/// \cond INTERNAL

namespace casadi{

  /** \brief Internal storage for DpleSolver related data

      @copydoc DPLE_doc
     \author Joris Gillis
      \date 2014
  */
  class CASADI_CONTROL_EXPORT DpleInternal : public FunctionInternal{
  public:
    /** \brief  Constructor
     *  \param[in] A  List of sparsities of A_i
     *  \param[in] V  List of sparsities of V_i
     */
    DpleInternal(const std::vector< Sparsity > & A, const std::vector< Sparsity > &V,
                 int nfwd=0, int nadj=0);

    /** \brief  Destructor */
    virtual ~DpleInternal()=0;

    /** \brief  Clone */
    virtual DpleInternal* clone() const=0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual DpleInternal* create(const std::vector< Sparsity > & A,
                                 const std::vector< Sparsity > &V) const = 0;

    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const{}

    /** \brief  evaluate */
    virtual void evaluate()=0;

    /** \brief  Initialize */
    virtual void init();

    /** \brief Generate a function that calculates \a nfwd forward derivatives
     and \a nadj adjoint derivatives
     */
    virtual Function getDerivative(int nfwd, int nadj)=0;

    /// List of sparsities of A_i
    std::vector< Sparsity > A_;

    /// List of sparsities of V_i
    std::vector< Sparsity > V_;

    /// Period
    int K_;

    /// Constant dimensions
    bool const_dim_;

    /// Assume positive definiteness of P_i
    bool pos_def_;

    /// Throw an error when system is unstable
    bool error_unstable_;

    /// Margin for instability detection
    double eps_unstable_;

    /// Number of forward derivatives
    int nfwd_;

    /// Number of adjoint derivatives
    int nadj_;

  };

} // namespace casadi
/// \endcond
#endif // DPLE_INTERNAL_HPP
