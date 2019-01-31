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

#ifndef CASADI_DM_HPP
#define CASADI_DM_HPP

#include "dm_fwd.hpp"
#include "matrix_decl.hpp"

namespace casadi {


  template<>
  DM DM::
  solve(const DM& A, const DM& b,
        const std::string& lsolver, const Dict& dict);

  template<>
  DM DM::
  inv(const DM& A,
        const std::string& lsolver, const Dict& dict);
  template<>
  DM DM::
  pinv(const DM& A, const std::string& lsolver,
       const Dict& dict);

  template<>
  DM DM::
  rand(const Sparsity& sp); // NOLINT(runtime/threadsafe_fn)

  template<>
  DM DM::
  expm(const DM& A);

  template<>
  DM DM::
  expm_const(const DM& A, const DM& t);

  template<> void DM::export_code(const std::string& lang,
       std::ostream &stream, const Dict& options) const;

  template<>
  Dict DM::info() const;

  template<>
  void DM::to_file(const std::string& filename, const Sparsity& sp,
    const double* nonzeros, const std::string& format);

  template<>
  DM DM::from_file(const std::string& filename, const std::string& format_hint);

#ifndef CASADI_DM_INSTANTIATOR_CPP
  extern template class Matrix<double>;
#endif // CASADI_DM_INSTANTIATOR_CPP

} // namespace casadi

#endif // CASADI_DM_HPP
