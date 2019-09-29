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

#define CASADI_IM_INSTANTIATOR_CPP
#include "matrix_impl.hpp"

using namespace std;

namespace casadi {

  bool CASADI_EXPORT is_slice(const IM& x, bool ind1) {
    return x.is_scalar() || (x.is_column() && x.is_dense() && is_slice(x.nonzeros(), ind1));
  }

  Slice CASADI_EXPORT to_slice(const IM& x, bool ind1) {
    return x.is_scalar() ? Slice(x.scalar(), ind1) : to_slice(x.nonzeros(), ind1);
  }

  template<>
  Dict CASADI_EXPORT IM::info() const {
    return {{"sparsity", sparsity().info()}, {"data", nonzeros()}};
  }
  template<>
  void CASADI_EXPORT IM::to_file(const std::string& filename,
      const Sparsity& sp, const casadi_int* nonzeros,
      const std::string& format_hint) {
    casadi_error("Not implemented");
  }

  // Instantiate templates
  template class casadi_limits<casadi_int>;
  template class CASADI_EXPORT Matrix<casadi_int>;


} // namespace casadi
