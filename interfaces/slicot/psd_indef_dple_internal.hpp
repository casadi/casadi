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

#ifndef PSD_INDEF_DPLE_INTERNAL_HPP
#define PSD_INDEF_DPLE_INTERNAL_HPP

#include "psd_indef_dple_solver.hpp"
#include "../../control/dple_internal.hpp"

namespace CasADi{

  /** \brief Internal storage for DpleSolver related data

      @copydoc DPLE_doc
     \author Joris gillis
      \date 2014
  */
  class PsdIndefDpleInternal : public DpleInternal{
  public:
    /** \brief  Constructor
     *  \param[in] A  List of sparsities of A_i 
     *  \param[in] V  List of sparsities of V_i 
     */
    PsdIndefDpleInternal(const std::vector< CRSSparsity > & A, const std::vector< CRSSparsity > &V, int nwfd=0, int nadj=0);
    
    /** \brief  Destructor */
    virtual ~PsdIndefDpleInternal();

    /** \brief  Clone */
    virtual PsdIndefDpleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
  
    /** \brief  Create a new solver */
    virtual PsdIndefDpleInternal* create(const std::vector< CRSSparsity > & A, const std::vector< CRSSparsity > &V) const{ return new PsdIndefDpleInternal(A,V);}
     
    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const{}

    /** \brief  evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    /// Generate a function that calculates nfwd forward derivatives and nadj adjoint derivatives
    virtual FX getDerivative(int nfwd, int nadj);
    
  private:
    /// Dimension of state-space
    int n_;
    
    /// Hessenberg-triangular data
    std::vector<double> T_;

    /// Schur form multiplier data
    std::vector<double> Z_;
    
    /// Schur form multiplier data
    std::vector<double> X_;
    
    /// Schur form multiplier data
    std::vector<double> dX_;
    
    // Schur form multiplier data
    std::vector<double> Xbar_;
    
    /// Transformed V data
    std::vector<double> VZ_;
    
    /// Temp data  nxn x K
    std::vector< Matrix<double> > nnKa_;
    
    /// Temp data  nxn x K
    std::vector< Matrix<double> > nnKb_;
  
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
    
  };
  
} // namespace CasADi

#endif // PSD_INDEF_DPLE_INTERNAL_HPP
