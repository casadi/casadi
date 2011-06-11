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

#ifndef NONZEROS_HPP
#define NONZEROS_HPP

namespace CasADi{


/** NonZeros class for Matrix 
  NonZeros is the return type for operator[] of the Matrix class, it allows access to the value as well as changing the parent object
  \author Joel Andersson 
  \date 2011
*/

/// submatrix - move to a new file
template<typename M>
class NonZeros : public M{
  public:
    /// Constructor
    NonZeros(M& mat, const std::vector<int>& kk) : M(mat.getNZ(kk)), mat_(mat), kk_(kk){}

    //@{
    /// Methods that modify a part of the parent obejct (A[k] = ?, A[k] += ?, etc.)
    const M& operator=(const NonZeros<M> &y);
    const M& operator=(const M &y);
    M operator+=(const M &y);
    M operator-=(const M &y);
    M operator*=(const M &y);
    M operator/=(const M &y);
    //@}

  private:
    /// A reference to the matrix that is allowed to be modified
    M& mat_;
    
    /// The element of the matrix that is allowed to be modified
    std::vector<int> kk_;
};

// Implementation
template<typename M>
const M& NonZeros<M>::operator=(const NonZeros<M> &y){ 
  mat_.setNZ(kk_,y); 
  return y;
}

// Implementation
template<typename M>
const M& NonZeros<M>::operator=(const M &y) { 
  mat_.setNZ(kk_,y); 
  return y;
}

template<typename M>
M NonZeros<M>::operator+=(const M &y){ 
  M s = *this+y;
  mat_.setNZ(kk_,s); 
  return s;
}

template<typename M>
M NonZeros<M>::operator-=(const M &y){ 
  M s = *this-y;
  mat_.setNZ(kk_,s); 
  return s;
}

template<typename M>
M NonZeros<M>::operator*=(const M &y){ 
   M s = *this*y;
   mat_.setNZ(kk_,s); 
   return mat_;
}

template<typename M>
M NonZeros<M>::operator/=(const M &y){ 
  M s = *this/y;
  mat_.setNZ(kk_,s); 
  return s;
}


} // namespace CasADi


#endif // NONZEROS_HPP
