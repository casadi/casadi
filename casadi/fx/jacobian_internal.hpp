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
    
    /// New constructor (not yet working)
    JacobianInternal(const FX& fcn, const std::vector<std::pair<int,int> >& jblocks);

    /// Clone
    virtual JacobianInternal* clone() const;

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
    
    /// Destructor
    virtual ~JacobianInternal();
      
    /// Evaluate the jacobian
    virtual void evaluate(int nfdir, int nadir);

    /// Initialize
    virtual void init();
  
    // Function to be differentiated
    FX fcn_;
  
    // Jacobian blocks requested
    std::vector<std::pair<int,int> > jblocks_;

    // Seeding matrices
    std::vector<CRSSparsity> D1_, D2_;

    // Number of simultaineous forward or adjoint directions of the function to be differentiated
    int nadir_fcn_;
    int nfdir_fcn_;

    // Transpose of the Jacobian sparsity
    CRSSparsity js_;
    CRSSparsity js_trans_;
    std::vector<int> js_trans_mapping_;
    
    // Output and input indices
    int iind_, oind_;
};

} // namespace CasADi


#endif // JACOBIAN_INTERNAL_HPP

