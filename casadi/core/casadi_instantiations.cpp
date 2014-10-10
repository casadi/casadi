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


#include "sx/sx_element.hpp"
#include "matrix/matrix.hpp"
#include "matrix/sparse_storage.hpp"
#include "mx/mx.hpp"

#include "casadi_limits.hpp"

#include "weak_ref.hpp"
#include <iostream>

#include "function/schemes_helpers.hpp"

#include "matrix/sparse_storage_impl.hpp"

using namespace std;
namespace casadi {

  INSTANTIATE_SUBMATRIX(Matrix<SXElement>)
  INSTANTIATE_NONZEROS(Matrix<SXElement>)

  template class GenericMatrix< Matrix<SXElement> >;
  template class Matrix< SXElement >;

  INSTANTIATE_SUBMATRIX(MX)
  INSTANTIATE_NONZEROS(MX)

  INSTANTIATE_SUBMATRIX(Matrix<int>)
  INSTANTIATE_SUBMATRIX(Matrix<double>)
  INSTANTIATE_NONZEROS(Matrix<int>)
  INSTANTIATE_NONZEROS(Matrix<double>)

  template class GenericMatrix< Matrix<double> >;
  template class GenericMatrix< Matrix<int> >;

  template class Matrix<double>;
  template class Matrix<int>;

  template<class T>
  const T casadi_limits<T>::zero = T(0);

  template<class T>
  const T casadi_limits<T>::one = 1;

  template<class T>
  const T casadi_limits<T>::two = 2;

  template<class T>
  const T casadi_limits<T>::minus_one = -1;

  template class casadi_limits<double>;
  template class casadi_limits<int>;

  template class SparseStorage<Sparsity>;
  template class SparseStorage<WeakRef>;


  INSTANTIATE_IOSCHEME_HELPERS(SX)
  INSTANTIATE_IOSCHEME_HELPERS(MX)
  INSTANTIATE_IOSCHEME_HELPERS(Sparsity)

} // namespace casadi

namespace std {
  #ifndef _MSC_VER
  template class std::numeric_limits<casadi::SXElement>;
  #endif

} // namespace std
