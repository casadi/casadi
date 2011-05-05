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

#ifndef CSPARSE_TOOLS_HPP
#define CSPARSE_TOOLS_HPP

#include "casadi/matrix/matrix.hpp"

// Wrapper around some useful functions in csparse

namespace CasADi{
  namespace Interfaces{
    
    /// BLT decomposition
    class BLT{
      public:
        /// Constructor
        BLT(CRSSparsity s, int seed=0);
        
        /// BLT transformation
        std::vector<int> rowperm; // row permutations
        std::vector<int> colperm; // column permutations
        std::vector<int> rowblock;  // block k is rows r[k] to r[k+1]-1
        std::vector<int> colblock;  // block k is cols s[k] to s[k+1]-1
        int nb;
        std::vector<int> coarse_rowblock;  // coarse row decomposition
        std::vector<int> coarse_colblock;  //coarse column decomposition
    };
    
  } // namespace Interfaces
} // namespace CasADi

#endif //CSPARSE_INTERNAL_HPP

