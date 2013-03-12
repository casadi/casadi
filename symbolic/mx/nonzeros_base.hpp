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

#ifndef NONZEROS_BASE_HPP
#define NONZEROS_BASE_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

namespace CasADi{
  /** \brief Abstract base class for GetNonzeros, AddNonzeros and SetNonzeros
      \author Joel Andersson
      \date 2013
  */
  class NonzerosBase : public MXNode{
  public:

    /// Constructor
    NonzerosBase(const std::vector<int>& nz) : nz_(nz){}

    /// Destructor
    virtual ~NonzerosBase() = 0;

    /// Clone function
    virtual NonzerosBase* clone() const = 0;
    
    /// Common implementation of AddNonzeros and SetNonzeros
    void evaluateMXBase(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given, bool add);

    /// Operation sequence
    std::vector<int> nz_;
  };

} // namespace CasADi

#endif // NONZEROS_BASE_HPP
