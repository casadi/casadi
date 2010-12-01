/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef ELEMENT_HPP
#define ELEMENT_HPP

namespace CasADi{


/** Element class for Matrix 
  Element is the return type for operator() of the Matrix class.
  From several alternative solutions, this was the only one which didn't cause ambigousity.
  Suggestions for improvement are welcome.
  \author Joel Andersson 
  \date 2010
*/

template<typename M, typename T>
class Element{
  public:
    /// Constructor
    Element(M& mat, int i, int j);
    
    //@{
    /// Methods that modify a part of the parent obejct (A[i] = ?, A[i] += ?, etc.)
    M& operator=(const T &y);
    M& operator+=(const T &y);
    M& operator-=(const T &y);
    M& operator*=(const T &y);
    M& operator/=(const T &y);
    //@}

    /// Get a reference to the element (? = A[i], ? += A[i], etc.)
    operator T();
  
  private:
    M& mat_;
    int i_, j_;
};


// Implementation

template<typename M, typename T>
M& Element<M,T>::operator=(const T &y){
  mat_.getElementRef(i_,j_) = y;
  return mat_;
}

template<typename M, typename T>
M& Element<M,T>::operator+=(const T &y){
  mat_.getElementRef(i_,j_) += y;
  return mat_;
}

template<typename M, typename T>
M& Element<M,T>::operator-=(const T &y){
  mat_.getElementRef(i_,j_) -= y;
  return mat_;
}

template<typename M, typename T>
M& Element<M,T>::operator*=(const T &y){
  mat_.getElementRef(i_,j_) *= y;
  return mat_;
}

template<typename M, typename T>
M& Element<M,T>::operator/=(const T &y){
  mat_.getElementRef(i_,j_) /= y;
  return mat_;
}

template<typename M, typename T>
Element<M,T>::operator T(){
  return mat_.getElement(i_,j_);
}

template<typename M, typename T>
Element<M,T>::Element(M& mat, int i, int j) : mat_(mat), i_(i), j_(j){
}



} // namespace CasADi


#endif // ELEMENT_HPP
