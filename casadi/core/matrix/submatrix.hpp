/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_SUBMATRIX_HPP
#define CASADI_SUBMATRIX_HPP

namespace casadi {


  /** SubMatrix class for Matrix
      SubMatrix is the return type for operator() of the Matrix class, it allows access to the value as well as changing the parent object
      \author Joel Andersson
      \date 2011
  */

  /// submatrix
  template<typename M, typename I, typename J>
  class CASADI_CORE_EXPORT SubMatrix : public M {
  public:
    /// Constructor
    SubMatrix(M& mat, const I& i, const J& j) : M(mat.sub(i, j)), mat_(mat), i_(i), j_(j) {}

    ///@{
    /// Methods that modify a part of the parent object (A(i, j) = ?, A(i, j) += ?, etc.)
    const M& operator=(const SubMatrix<M, I, J> &y);
    const M& operator=(const M &y);
    M operator+=(const M &y);
    M operator-=(const M &y);
    M operator*=(const M &y);
    M operator/=(const M &y);
    ///@}

  private:
    /// A reference to the matrix that is allowed to be modified
    M& mat_;

    /// The element of the matrix that is allowed to be modified
    I i_;
    J j_;
  };

#ifdef casadi_core_implementation
  // Implementation
  template<typename M, typename I, typename J>
  const M& SubMatrix<M, I, J>::operator=(const SubMatrix<M, I, J> &y) {
    mat_.setSub(y, i_, j_);
    return y;
  }

  // Implementation
  template<typename M, typename I, typename J>
  const M& SubMatrix<M, I, J>::operator=(const M &y) {
    mat_.setSub(y, i_, j_);
    return y;
  }

  template<typename M, typename I, typename J>
  M SubMatrix<M, I, J>::operator+=(const M &y) {
    M s = *this+y;
    mat_.setSub(s, i_, j_);
    return s;
  }

  template<typename M, typename I, typename J>
  M SubMatrix<M, I, J>::operator-=(const M &y) {
    M s = *this-y;
    mat_.setSub(s, i_, j_);
    return s;
  }

  template<typename M, typename I, typename J>
  M SubMatrix<M, I, J>::operator*=(const M &y) {
    M s = *this*y;
    mat_.setSub(s, i_, j_);
    return s;
  }

  template<typename M, typename I, typename J>
  M SubMatrix<M, I, J>::operator/=(const M &y) {
    M s = *this/y;
    mat_.setSub(s, i_, j_);
    return s;
  }
#endif

#define INSTANTIATE_SUBMATRIX(Mt)\
template class SubMatrix< Mt , int, int >;\
template class SubMatrix< Mt , int, std::vector<int> >;\
template class SubMatrix< Mt , int, Slice >;\
template class SubMatrix< Mt , std::vector<int> , int >;\
template class SubMatrix< Mt , std::vector<int> , std::vector<int> >;\
template class SubMatrix< Mt , std::vector<int> , Matrix<int> >;\
template class SubMatrix< Mt , std::vector<int> , Slice >;\
template class SubMatrix< Mt , int , Matrix<int> >;\
template class SubMatrix< Mt , Slice, int >;\
template class SubMatrix< Mt , Slice, Slice >;\
template class SubMatrix< Mt , Matrix<int> , int >;\
template class SubMatrix< Mt , Matrix<int> , std::vector<int> >;\
template class SubMatrix< Mt , Matrix<int> , Matrix<int> >;\
template class SubMatrix< Mt , Sparsity , int >;

} // namespace casadi


#endif // CASADI_SUBMATRIX_HPP
