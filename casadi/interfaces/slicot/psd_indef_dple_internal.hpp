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

#ifndef CASADI_PSD_INDEF_DPLE_INTERNAL_HPP
#define CASADI_PSD_INDEF_DPLE_INTERNAL_HPP

#include "../../control/dple_internal.hpp"
#include <casadi/interfaces/slicot/casadi_dplesolver_slicot_export.h>

/// \cond INTERNAL
namespace casadi {

  /** \brief An efficient solver for Discrete Periodic Lyapunov Equations using SLICOT
   *
   * @copydoc DPLE_doc
   
       Uses Periodic Schur Decomposition ('psd') and does not assume positive definiteness.
       Based on Periodic Lyapunov equations: some applications and new algorithms.
       Int. J. Control, vol. 67, pp. 69-87, 1997.

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLESOLVER_SLICOT_EXPORT PsdIndefDpleInternal : public DpleInternal {
  public:
    /** \brief  Constructor
     *  \param[in] A  List of sparsities of A_i
     *  \param[in] V  List of sparsities of V_i
     */
    PsdIndefDpleInternal(const std::vector< Sparsity > & A, const std::vector< Sparsity > &V,
                         int nrhs=1, bool transp=false);

    /** \brief  Destructor */
    virtual ~PsdIndefDpleInternal();

    /** \brief  Clone */
    virtual PsdIndefDpleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual PsdIndefDpleInternal* create(const std::vector< Sparsity >& A,
                                         const std::vector< Sparsity >& V) const
    { return new PsdIndefDpleInternal(A, V); }

    /** \brief  Create a new DPLE Solver */
    static DpleInternal* creator(const std::vector< Sparsity >& A,
                                 const std::vector< Sparsity >& V)
    { return new PsdIndefDpleInternal(A, V);}

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

  private:
    /// Dimension of state-space
    int n_;

    /// Hessenberg-triangular data
    std::vector<double> T_;

    /// Schur form multiplier data
    std::vector<double> Z_;

    /// Schur form multiplier data
    std::vector<double> X_;

    // Schur form multiplier data
    std::vector<double> Xbar_;

    /// Transformed V data
    std::vector<double> VZ_;

    /// Temp data  (n x n) x K
    std::vector< Matrix<double> > nnKa_;

    /// Temp data  (n x n) x K
    std::vector< Matrix<double> > nnKb_;

    /// Real parts of eigenvalues
    std::vector< double > eig_real_;

    /// Imaginary parts of eigenvalues
    std::vector< double > eig_imag_;

    /// Solvers for low-order Discrete Periodic Sylvester Equations
    std::vector< std::vector< LinearSolver> > dpse_solvers_;

    std::vector<int> partition_;

    int partindex(int i, int j, int k, int r, int c);

    /// Temp data  F
    std::vector< double > F_;

    /// Temp data  FF
    std::vector< double > FF_;

    /// Work vector for periodic Schur form
    std::vector< double > dwork_;

    /// Numerical zero, used in periodic Schur form
    double psd_num_zero_;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_PSD_INDEF_DPLE_INTERNAL_HPP
