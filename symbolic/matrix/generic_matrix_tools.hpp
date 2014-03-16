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
  template<typename MatType>
  MatType linspace(const GenericMatrix<MatType> &a, const GenericMatrix<MatType> &b, int nsteps);

  /** \brief Matlab's cross command
   */
  template<typename MatType>
  MatType cross(const GenericMatrix<MatType> &a, const GenericMatrix<MatType> &b, int dim = -1);

  /** \brief Convert a lower triangular matrix to a symmetric one
   */
  template<typename MatType>
  MatType tril2symm(const GenericMatrix<MatType> &a);

  /** \brief Convert a upper triangular matrix to a symmetric one
   */
  template<typename MatType>
  MatType triu2symm(const GenericMatrix<MatType> &a);

  /** \brief Get the upper triangular part of a matrix
   */
  template<typename MatType>
  MatType triu(const GenericMatrix<MatType> &a);

  /** \brief Get the lower triangular part of a matrix
   */
  template<typename MatType>
  MatType tril(const GenericMatrix<MatType> &a);

  /** \brief Check if two expressions are equal, assuming that they are comparible */
  template<typename MatType>
  bool isEqual(const GenericMatrix<MatType>& x, const GenericMatrix<MatType>& y){ 
    return static_cast<const MatType&>(x).isEqual(static_cast<const MatType&>(y));
  }

#ifndef SWIG
  template<typename MatType>
  MatType linspace(const GenericMatrix<MatType> &a_, const GenericMatrix<MatType> &b_, int nsteps){
    const MatType& a = static_cast<const MatType&>(a_);
    const MatType& b = static_cast<const MatType&>(b_);
    std::vector<MatType> ret(nsteps);
    ret[0] = a;
    MatType step = (b-a)/(nsteps-1);

    for(int i=1; i<nsteps-1; ++i)
      ret[i] = ret[i-1] + step;
  
    ret[nsteps-1] = b;
    return vertcat(ret);
  }
#endif // SWIG

#ifndef SWIG
  template<typename MatType>
  MatType cross(const GenericMatrix<MatType> &a, const GenericMatrix<MatType> &b, int dim) {
    casadi_assert_message(a.size1()==b.size1() && a.size2()==b.size2(),"cross(a,b): Inconsistent dimensions. Dimension of a (" << a.dimString() << " ) must equal that of b (" << b.dimString() << ").");
  
    casadi_assert_message(a.size1()==3 || a.size2()==3,"cross(a,b): One of the dimensions of a should have length 3, but got " << a.dimString() << ".");
    casadi_assert_message(dim==-1 || dim==1 || dim==2,"cross(a,b,dim): Dim must be 1, 2 or -1 (automatic).");
  
  
    std::vector<MatType> ret(3);
  
  
    bool t = a.size1()==3;
  
    if (dim==1) t = true;
    if (dim==2) t = false;
  
    MatType a1 = t ? a(0,ALL) : a(ALL,0);
    MatType a2 = t ? a(1,ALL) : a(ALL,1);
    MatType a3 = t ? a(2,ALL) : a(ALL,2);

    MatType b1 = t ? b(0,ALL) : b(ALL,0);
    MatType b2 = t ? b(1,ALL) : b(ALL,1);
    MatType b3 = t ? b(2,ALL) : b(ALL,2);
    
    ret[0] = a2*b3-a3*b2;
    ret[1] = a3*b1-a1*b3;
    ret[2] = a1*b2-a2*b1;
    
    return t ? vertcat(ret) : horzcat(ret);
  
  }

  template<typename MatType> 
  MatType tril2symm(const GenericMatrix<MatType> &a_) {
    const MatType& a = static_cast<const MatType&>(a_);
    casadi_assert_message(a.isSquare(),"Shape error in tril2symm. Expecting square shape but got " << a.dimString());
    casadi_assert_message(a.sizeU()-a.sizeD()==0,"Sparsity error in tril2symm. Found above-diagonal entries in argument: " << a.dimString());
    return a +  a.T() - diag(diag(a));
  }


  template<typename MatType> 
  MatType triu2symm(const GenericMatrix<MatType> &a_) {
    const MatType& a = static_cast<const MatType&>(a_);
    casadi_assert_message(a.isSquare(),"Shape error in triu2symm. Expecting square shape but got " << a.dimString());
    casadi_assert_message(a.sizeL()-a.sizeD()==0,"Sparsity error in triu2symm. Found below-diagonal entries in argument: " << a.dimString());
    return a + a.T() - diag(diag(a));
  }

  template<typename MatType>
  MatType triu(const GenericMatrix<MatType> &a_) {
    const MatType& a = static_cast<const MatType&>(a_);
    return a.setSparse(a.sparsity().getTriu());
  }

  template<typename MatType>
  MatType tril(const GenericMatrix<MatType> &a_) {
    const MatType& a = static_cast<const MatType&>(a_);
    return a.setSparse(a.sparsity().getTril());
  }
#endif // SWIG


} // namespace CasADi

#ifdef SWIG

// map the template name to the instantiated name
#define GMTT_INST(MatType,function_name) \
%template(function_name) CasADi::function_name< MatType >;

// Define template instanciations
#define GENERIC_MATRIX_TOOLS_TEMPLATES(MatType) \
GMTT_INST(MatType,cross) \
GMTT_INST(MatType,tril2symm) \
GMTT_INST(MatType,triu2symm) \
GMTT_INST(MatType,triu) \
GMTT_INST(MatType,tril) \
GMTT_INST(MatType,isEqual)

// Define template instanciations
#define GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(MatType) \
GMTT_INST(MatType,linspace)


#endif //SWIG



#endif // GENERIC_MATRIX_TOOLS_HPP

