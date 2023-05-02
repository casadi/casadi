/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

#ifndef CASADI_IM_HPP
#define CASADI_IM_HPP

#include "im_fwd.hpp"
#include "matrix_decl.hpp"

namespace casadi {

  /// Is the IM a Slice
  bool CASADI_EXPORT is_slice(const IM& x, bool ind1=false);

  ///  Convert IM to Slice
  Slice CASADI_EXPORT to_slice(const IM& x, bool ind1=false);

  template<>
  Dict CASADI_EXPORT IM::info() const;
  template<>
  void CASADI_EXPORT IM::to_file(const std::string& filename,
    const Sparsity& sp, const casadi_int* nonzeros,
    const std::string& format_hint);

#ifndef CASADI_IM_INSTANTIATOR_CPP
  extern template class Matrix<casadi_int>;
#endif // CASADI_IM_INSTANTIATOR_CPP

} // namespace casadi

#endif // CASADI_IM_HPP
