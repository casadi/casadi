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

#ifndef SLICE_HPP
#define SLICE_HPP

#include <vector>
#include "../casadi_exception.hpp"
#include <limits>

namespace CasADi{
  
  /// Dummy class denoting all rows/columns
  class Slice{
    public:
      /// Constructor
      Slice(int start__=0, int stop__=std::numeric_limits<int>::max(), int step__=1);
      
      /// Get a vector of indices
      std::vector<int> getAll(int len) const;
      
      /// Data members (all public)
      int start;
      int stop;
      int step;
  };
  static Slice ALL;
  
    /// Dummy class denoting all rows/columns
  class IndexList{
    private:

    public:
      enum Type {NILL, INT, SLICE, IVECTOR};
      /// Constructor
      IndexList();
      IndexList(int i);
      IndexList(const std::vector<int> &i);
      IndexList(const Slice &i);
      
      /// Get a vector of indices
      std::vector<int> getAll(int len) const;
      
      /// Data members (all public)
      Slice slice;
      int i;
      std::vector<int> iv;
      Type type;
  };
  
} // namespace CasADi


#ifdef SWIG
%template(Pair_Slice_Int) std::pair<CasADi::Slice,int>;
%template(Pair_Int_Slice) std::pair<int,CasADi::Slice>;
%template(Pair_Slice_Slice) std::pair<CasADi::Slice,CasADi::Slice>;
#endif // SWIG

#endif // SLICE_HPP

