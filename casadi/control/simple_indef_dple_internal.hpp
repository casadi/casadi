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

#ifndef SIMPLE_INDEF_DPLE_INTERNAL_HPP
#define SIMPLE_INDEF_DPLE_INTERNAL_HPP

#include "simple_indef_dple_solver.hpp"
#include "dple_internal.hpp"

/// \cond INTERNAL
namespace casadi{

  /** \brief Internal storage for DpleSolver related data

      @copydoc DPLE_doc
     \author Joris gillis
      \date 2014
  */
  class CASADI_CONTROL_EXPORT SimpleIndefDpleInternal : public DpleInternal{
  public:
    /** \brief  Constructor
     *  \param[in] A  List of sparsities of A_i
     *  \param[in] V  List of sparsities of V_i
     */
    SimpleIndefDpleInternal(const std::vector< Sparsity > & A, const std::vector< Sparsity > &V);

    /** \brief  Destructor */
    virtual ~SimpleIndefDpleInternal();

    /** \brief  Clone */
    virtual SimpleIndefDpleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual SimpleIndefDpleInternal* create(const std::vector< Sparsity > & A,
                                            const std::vector< Sparsity > &V) const {
        return new SimpleIndefDpleInternal(A,V);}

    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const{}

    /** \brief  evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    /// Generate a function that calculates nfwd forward derivatives and nadj adjoint derivatives
    virtual Function getDerivative(int nfwd, int nadj);

  private:
    /// Main implementation as MXFunction
    Function f_;

    /// State space dimension
    int n_;

  };

} // namespace casadi
/// \endcond
#endif // SIMPLE_INDEF_DPLE_INTERNAL_HPP
