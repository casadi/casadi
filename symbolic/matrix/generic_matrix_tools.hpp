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

#ifndef GENERIC_MATRIX_TOOLS_HPP
#define GENERIC_MATRIX_TOOLS_HPP

#include "slice.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "sparsity.hpp"
#include "../casadi_math.hpp"
#include "../casadi_exception.hpp"

namespace CasADi{


/** \brief Matlab's linspace command
*/
template<typename T>
T linspace(const GenericMatrix<T> &a, const GenericMatrix<T> &b, int nsteps);

/** \brief Matlab's cross command
*/
template<typename T>
T cross(const GenericMatrix<T> &a, const GenericMatrix<T> &b, int dim = -1);

/** \brief Convert a lower triangular matrix to a symmetric one
*/
template<typename T>
T tril2symm(const GenericMatrix<T> &a);

/** \brief Convert a upper triangular matrix to a symmetric one
*/
template<typename T>
T triu2symm(const GenericMatrix<T> &a);

#ifndef SWIG
template<typename T>
T linspace(const GenericMatrix<T> &a_, const GenericMatrix<T> &b_, int nsteps){
  const T& a = static_cast<const T&>(a_);
  const T& b = static_cast<const T&>(b_);
  std::vector<T> ret(nsteps);
  ret[0] = a;
  T step = (b-a)/(nsteps-1);

  for(int i=1; i<nsteps-1; ++i)
    ret[i] = ret[i-1] + step;
  
  ret[nsteps-1] = b;
  return vertcat(ret);
}
#endif // SWIG

#ifndef SWIG
template<typename T>
T cross(const GenericMatrix<T> &a, const GenericMatrix<T> &b, int dim) {
  casadi_assert_message(a.size1()==b.size1() && a.size2()==b.size2(),"cross(a,b): Inconsistent dimensions. Dimension of a (" << a.dimString() << " ) must equal that of b (" << b.dimString() << ").");
  
  casadi_assert_message(a.size1()==3 || a.size2()==3,"cross(a,b): One of the dimensions of a should have length 3, but got " << a.dimString() << ".");
  casadi_assert_message(dim==-1 || dim==1 || dim==2,"cross(a,b,dim): Dim must be 1, 2 or -1 (automatic).");
  
  
  std::vector<T> ret(3);
  
  
  bool t = a.size1()==3;
  
  if (dim==1) t = true;
  if (dim==2) t = false;
  
  T a1 = t ? a(0,ALL) : a(ALL,0);
  T a2 = t ? a(1,ALL) : a(ALL,1);
  T a3 = t ? a(2,ALL) : a(ALL,2);

  T b1 = t ? b(0,ALL) : b(ALL,0);
  T b2 = t ? b(1,ALL) : b(ALL,1);
  T b3 = t ? b(2,ALL) : b(ALL,2);
    
  ret[0] = a2*b3-a3*b2;
  ret[1] = a3*b1-a1*b3;
  ret[2] = a1*b2-a2*b1;
    
  return t ? vertcat(ret) : horzcat(ret);
  
}

template<typename T> 
T tril2symm(const GenericMatrix<T> &a_) {
  const T& a = static_cast<const T&>(a_);
  casadi_assert_message(a.isSquare(),"Shape error in tril2symm. Expecting square shape but got " << a.dimString());
  casadi_assert_message(a.sizeU()-a.sizeD()==0,"Sparsity error in tril2symm. Found above-diagonal entries in argument: " << a.dimString());
  return a + trans(a) - diag(diag(a));
}


template<typename T> 
T triu2symm(const GenericMatrix<T> &a_) {
  const T& a = static_cast<const T&>(a_);
  casadi_assert_message(a.isSquare(),"Shape error in triu2symm. Expecting square shape but got " << a.dimString());
  casadi_assert_message(a.sizeL()-a.sizeD()==0,"Sparsity error in triu2symm. Found below-diagonal entries in argument: " << a.dimString());
  return a + trans(a) - diag(diag(a));
}
#endif // SWIG


} // namespace CasADi

#ifdef SWIG

// map the template name to the instantiated name
#define GMTT_INST(T,function_name) \
%template(function_name) CasADi::function_name< T >;

// Define template instanciations
#define GENERIC_MATRIX_TOOLS_TEMPLATES(T) \
GMTT_INST(T,cross) \
GMTT_INST(T,tril2symm) \
GMTT_INST(T,triu2symm)

// Define template instanciations
#define GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(T) \
GMTT_INST(T,linspace)


#endif //SWIG



#endif // GENERIC_MATRIX_TOOLS_HPP

