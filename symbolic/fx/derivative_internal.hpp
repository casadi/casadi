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

#ifndef DERIVATIVE_INTERNAL_HPP
#define DERIVATIVE_INTERNAL_HPP

#include <vector>
#include "derivative.hpp"
#include "fx_internal.hpp"

namespace CasADi{
 
  /** \brief  Internal node class for Derivative
  \author Joel Andersson 
  \date 2012
*/
class DerivativeInternal : public FXInternal{
  friend class Derivative;
  public:
    
    /// New constructor
    DerivativeInternal(const FX& fcn, int nfwd, int nadj);

    /// Clone
    virtual DerivativeInternal* clone() const;

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
    
    /// Destructor
    virtual ~DerivativeInternal();
      
    /// Evaluate the jacobian
    virtual void evaluate(int nfdir, int nadir);

    /// Initialize
    virtual void init();
  
    // Function to be differentiated
    FX fcn_;
  
    // Number of directional derivatives
    int nfwd_, nadj_;
};

} // namespace CasADi


#endif // DERIVATIVE_INTERNAL_HPP

