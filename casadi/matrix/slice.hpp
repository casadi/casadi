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
  
  /** Class representing a slice
  */
  class Slice{
    public:
      /// Defailt constructor - all elements
      Slice();
      
      /// A single element
      Slice(int i);
      
      /// A slice
      Slice(int start, int stop, int step=1);
      
      /// Get a vector of indices
      std::vector<int> getAll(int len) const;
      
      /// Data members (all public)
      int start_;
      int stop_;
      int step_;
  };
  static Slice ALL;
  
  /** Class representing a set of indices of arbitrary order */
  class IndexSet{
    public:
    
      /// A single element
      IndexSet(int i);
  
      /// A set of indices
      IndexSet(const std::vector<int>& v);
      
      /// Get a vector of indices
      const std::vector<int>& getAll(int len) const;
      
      /// Data members (all public)
      std::vector<int> v_;
  };
  
  
   /**  Class representing a non-regular (and thus non-slice) index list 
   */
  class IndexList{
    private:

    public:
      enum Type {NILL, INT, SLICE, IVECTOR};
      /// Constructor
      IndexList();
      explicit IndexList(int i);
      explicit IndexList(const std::vector<int> &i);
      explicit IndexList(const Slice &i);
      
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

