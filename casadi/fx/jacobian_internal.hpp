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

#ifndef JACOBIAN_INTERNAL_HPP
#define JACOBIAN_INTERNAL_HPP

#include <vector>
#include "jacobian.hpp"
#include "fx_internal.hpp"

namespace CasADi{
 
  /** \brief  Internal node class for Jacobian
  \author Joel Andersson 
  \date 2010
*/
class JacobianInternal : public FXInternal{
  friend class Jacobian;
  public:
    
    /// Constructor
    JacobianInternal(const FX& fcn, int iind, int oind);

    /// Clone
    virtual JacobianInternal* clone() const;

    /// Destructor
    virtual ~JacobianInternal();
      
    /// Evaluate the jacobian
    virtual void evaluate(int nfdir, int nadir);

    /// Initialize
    virtual void init();
  
      
      
      /// Compress the number of evaluations using curtis-powel reed seeding
/*      void compress();*/
      
      bool use_fd_, use_ad_fwd_;
  
      std::vector<double> epsilon_; // perturbations
  
      // Dimensions
      int n_,m_;
  
      // Function to be differentiated
      FX fcn_;
  
      // Output and input indices
      int iind_, oind_;

      // Number of forward directions of the function to be differentiated
      int nadir_fcn_;
      int nfdir_fcn_;
      
      // Seeding, contains which jacobian row/column is calculated by which component seeding vector
      std::vector<std::vector<int> > seeding_;
};



} // namespace CasADi


#endif // JACOBIAN_INTERNAL_HPP

